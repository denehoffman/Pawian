//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <map>
#include <iterator>
#include <memory>
#include <memory>

#include "TROOT.h"

#include "epemUtils/epemParser.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/StreamFitParmsBase.hh"
#include "PwaUtils/PwaFcnBase.hh"
#include "PwaUtils/PwaCovMatrix.hh"
#include "PwaUtils/WaveContribution.hh"
#include "PwaUtils/AppBase.hh"

#include "Utils/PawianCollectionUtils.hh"
#include "Utils/ErrLogUtils.hh"
#include "epemUtils/epemEnv.hh"
#include "epemUtils/epemReaction.hh"
#include "epemUtils/epemBaseLh.hh"

#include "Event/EventReaderDefault.hh"
#include "epemUtils/epemHist.hh"

#include "Event/Event.hh"
#include "Event/EventList.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnScan.h"


#ifdef _OPENMP
#include <omp.h>
#endif


int main(int __argc,char *__argv[]){
  clock_t start, end;
  start= clock();

  // Disable output buffering
  setvbuf(stdout, NULL, _IONBF, 0);

  // Parse the command line
  epemParser* theAppParams=new epemParser(__argc, __argv);

  // Set the desired error logging mode
  setErrLogMode(theAppParams->getErrLogMode());
  
#ifdef _OPENMP
  const int noOfThreads=theAppParams->noOfThreads();
  omp_set_num_threads(noOfThreads);
#endif


  epemEnv::instance()->setup(theAppParams);

  std::shared_ptr<epemReaction> theEpEmReaction=epemEnv::instance()->reaction();

  theEpEmReaction->print(std::cout);

  std::string mode=theAppParams->mode();

  std::shared_ptr<FitParamsBase> theFitParamBase=std::shared_ptr<FitParamsBase>(new FitParamsBase());


  std::string prodFormalism=theAppParams->productionFormalism();
  std::shared_ptr<AbsLh> theLhPtr;
  theLhPtr=std::shared_ptr<AbsLh>(new epemBaseLh());

  AppBase theAppBase(epemEnv::instance(), theLhPtr, theFitParamBase);

  if (mode=="dumpDefaultParams"){
    theAppBase.dumpDefaultParams();
     return 0;
  }


  std::string paramStreamerPath=theAppParams->fitParamFile();
  std::string outputFileNameSuffix= epemEnv::instance()->outputFileNameSuffix();
  StreamFitParmsBase theParamStreamer(paramStreamerPath, theLhPtr);
  fitParams theStartparams=theParamStreamer.getFitParamVal();
  fitParams theErrorparams=theParamStreamer.getFitParamErr();

  if (mode=="gen"){
    theAppBase.generate(theStartparams);
    return 0;
  }


  const std::string datFile=theAppParams->dataFile();
  const std::string mcFile=theAppParams->mcFile();
  Info << "data file: " << datFile ;  // << endmsg;
  Info << "mc file: " << mcFile ;  // << endmsg;
  
  std::vector<std::string> dataFileNames;
  dataFileNames.push_back(datFile);

  std::vector<std::string> mcFileNames;
  mcFileNames.push_back(mcFile);  

  bool withEvtWeight=theAppParams->useEvtWeight();
  Info << "EvtWeight: " << withEvtWeight << endmsg;  

  EventList eventsData;
  theAppBase.readEvents(eventsData, dataFileNames, withEvtWeight);
  
  EventList mcData;
  theAppBase.readEvents(mcData, mcFileNames, withEvtWeight);


  std::shared_ptr<EvtDataBaseList> eventListPtr(new EvtDataBaseList(epemEnv::instance()));
  //  eventListPtr->ratioMcToData(theAppParams->ratioMcToData());
  eventListPtr->read(eventsData, mcData);

  theLhPtr->setDataVec(eventListPtr->getDataVecs());
  theLhPtr->setMcVec(eventListPtr->getMcVecs()); 

  PwaFcnBase theFcn(theLhPtr, theFitParamBase, outputFileNameSuffix);
  MnUserParameters upar;
  theFitParamBase->setMnUsrParams(upar, theStartparams, theErrorparams);
  
  std::cout << "\n\n**************** Minuit Fit parameter **************************" << std::endl;
  for (int i=0; i<int(upar.Params().size()); ++i){
    std::cout << upar.Name(i) << "\t" << upar.Value(i) << "\t" << upar.Error(i) << std::endl;
  }

  const std::vector<std::string> fixedParams=theAppParams->fixedParams();  
  const unsigned int noOfFreeFitParams = upar.Params().size()-fixedParams.size();

  if (mode=="qaMode"){
    double evtWeightSumData = eventListPtr->NoOfWeightedDataEvts();
    theAppBase.qaMode(theStartparams, evtWeightSumData, noOfFreeFitParams );
    epemHist theHist(theLhPtr, theStartparams);
    end= clock();
    double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
    Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;
    
    return 1;    
  }
  

  if (mode=="pwa"){
    bool cacheAmps = theAppParams->cacheAmps();
    Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
    if (cacheAmps) theLhPtr->cacheAmplitudes();
    std::vector<std::string>::const_iterator itFix;
    for (itFix=fixedParams.begin(); itFix!=fixedParams.end(); ++itFix){
      upar.Fix( (*itFix) );
    }

    MnMigrad migrad(theFcn, upar);
    Info <<"start migrad "<< endmsg;
    FunctionMinimum min = migrad();
    
    if(!min.IsValid()) {
      //try with higher strategy
      Info <<"FM is invalid, try with strategy = 2."<< endmsg;
      MnMigrad migrad2(theFcn, min.UserState(), MnStrategy(2));
      min = migrad2();
    }
    
    MnUserParameters finalUsrParameters=min.UserParameters();
    const std::vector<double> finalParamVec=finalUsrParameters.Params();
    fitParams finalFitParams=theStartparams;
    theFitParamBase->getFitParamVal(finalParamVec, finalFitParams);

    theFitParamBase->printParams(finalFitParams);
    double theLh=theLhPtr->calcLogLh(finalFitParams);
    Info <<"theLh = "<< theLh << endmsg;
    
    
    const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
    fitParams finalFitErrs=theErrorparams;
    theFitParamBase->getFitParamVal(finalParamErrorVec, finalFitErrs);
    
    std::ostringstream finalResultname;
    finalResultname << "finalResult" << outputFileNameSuffix << ".dat";

    std::ofstream theStream ( finalResultname.str().c_str() );
    theFitParamBase->dumpParams(theStream, finalFitParams, finalFitErrs);
    
    MnUserCovariance theCovMatrix = min.UserCovariance();
    std::cout  << min << std::endl;

    std::ostringstream serializationFileName;
    serializationFileName << "serializedOutput" << outputFileNameSuffix << ".dat";
    std::ofstream serializationStream(serializationFileName.str().c_str());
    boost::archive::text_oarchive boostOutputArchive(serializationStream);
    
    if(min.HasValidCovariance()){
      PwaCovMatrix thePwaCovMatrix(theCovMatrix, finalUsrParameters, finalFitParams);
      boostOutputArchive << thePwaCovMatrix;
    }
    
    return 0;
 }

}     
 
  
