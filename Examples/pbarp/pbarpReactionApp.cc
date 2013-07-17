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

#include "TROOT.h"

#include "pbarpUtils/pbarpParser.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include "pbarpUtils/pbarpStatesLS.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/StreamFitParmsBase.hh"
#include "PwaUtils/PwaFcnBase.hh"
#include "PwaUtils/PwaFcnServer.hh"
#include "PwaUtils/PwaCovMatrix.hh"
#include "PwaUtils/WaveContribution.hh"
//#include "PwaUtils/PwaGen.hh"
#include "PwaUtils/AppBase.hh"

#include "Utils/PawianCollectionUtils.hh"
#include "Utils/ErrLogUtils.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "pbarpUtils/pbarpReaction.hh"
#include "pbarpUtils/pbarpBaseLh.hh"
#include "pbarpUtils/pbarpHeliLh.hh"
#include "pbarpUtils/pbarpCanoLh.hh"
#include "pbarpUtils/pbarpTensorLh.hh"

#include "Event/EventReaderDefault.hh"

#include "PwaUtils/EvtDataBaseList.hh"
#include "pbarpUtils/pbarpHist.hh"
#include "pbarpUtils/spinDensityHist.hh"
#include "Event/Event.hh"
#include "Event/EventList.hh"
#include "PwaUtils/NetworkServer.hh"
#include "PwaUtils/NetworkClient.hh"
#include "PwaUtils/WelcomeScreen.hh"

#include "PwaUtils/EvoMinimizer.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnScan.h"



int main(int __argc,char *__argv[]){
  clock_t start, end;
  start= clock();

  Info << welcomeScreen << endmsg;

  // Disable output buffering
  setvbuf(stdout, NULL, _IONBF, 0);

  // Parse the command line
  pbarpParser* theAppParams=new pbarpParser(__argc, __argv);

  // Set the desired error logging mode
  setErrLogMode(theAppParams->getErrLogMode());
  

  pbarpEnv::instance()->setup(theAppParams);

  std::shared_ptr<pbarpReaction> thepbarpReaction=pbarpEnv::instance()->reaction();

  thepbarpReaction->print(std::cout);


  std::string mode=theAppParams->mode();

  std::shared_ptr<FitParamsBase> theFitParamBase=std::shared_ptr<FitParamsBase>(new FitParamsBase());


  std::string prodFormalism=theAppParams->productionFormalism();
  std::shared_ptr<AbsLh> theLhPtr;
  if(prodFormalism=="Cano") theLhPtr=std::shared_ptr<AbsLh>(new pbarpCanoLh());
  else if(prodFormalism=="Heli") theLhPtr=std::shared_ptr<AbsLh>(new pbarpHeliLh());
  else if(prodFormalism=="Tensor") theLhPtr=std::shared_ptr<AbsLh>(new pbarpTensorLh());
  else {
    Alert << "prodFormalism\t" << prodFormalism << "\tdoesn't exist!!!" << endmsg;
    exit(1);
  }

  AppBase theAppBase(pbarpEnv::instance(), theLhPtr, theFitParamBase);

  if (mode=="dumpDefaultParams"){
    theAppBase.dumpDefaultParams();
    return 1;
  }


  std::string paramStreamerPath=theAppParams->fitParamFile();
  std::string outputFileNameSuffix= pbarpEnv::instance()->outputFileNameSuffix();
  StreamFitParmsBase theParamStreamer(paramStreamerPath, theLhPtr);
  fitParams theStartparams=theParamStreamer.getFitParamVal();
  fitParams theErrorparams=theParamStreamer.getFitParamErr();

  if (mode=="gen"){
    theAppBase.generate(theStartparams);
    return 1;
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


  MnUserParameters upar;
  theFitParamBase->setMnUsrParams(upar, theStartparams, theErrorparams);
  
  std::cout << "\n\n**************** Minuit Fit parameter **************************" << std::endl;
  for (int i=0; i<int(upar.Params().size()); ++i){
    std::cout << upar.Name(i) << "\t" << upar.Value(i) << "\t" << upar.Error(i) << std::endl;
  }

  const std::vector<std::string> fixedParams=theAppParams->fixedParams();  
  const unsigned int noOfFreeFitParams = upar.Params().size()-fixedParams.size();



  if(mode == "client"){

  bool cacheAmps = theAppParams->cacheAmps();
  Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
  if (cacheAmps) theLhPtr->cacheAmplitudes();
  
  std::ostringstream portStringStream;
  portStringStream << theAppParams->serverPort();
  
  NetworkClient theClient(theAppParams->serverAddress(), portStringStream.str());
  if(!theClient.Login())
    return 0;
  
  
  EventList eventsDataClient;  
  theAppBase.readEvents(eventsDataClient, dataFileNames, withEvtWeight, theClient.GetEventLimits()[0], theClient.GetEventLimits()[1]);  
  
  EventList mcDataClient; 
  theAppBase.readEvents(mcDataClient, mcFileNames, false, theClient.GetEventLimits()[2], theClient.GetEventLimits()[3]);  
 
  std::shared_ptr<EvtDataBaseList> pbarpEventListPtr(new EvtDataBaseList(pbarpEnv::instance()));
  pbarpEventListPtr->read(eventsDataClient, mcDataClient);

  theLhPtr->setDataVec(pbarpEventListPtr->getDataVecs());
  theLhPtr->setMcVec(pbarpEventListPtr->getMcVecs());

  theAppBase.calcAndSendClientLh(theClient, theStartparams);
  
  return 1;
 }

  if(mode == "spinDensity"){

  bool cacheAmps = theAppParams->cacheAmps();
  Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
  if (cacheAmps) theLhPtr->cacheAmplitudes();

  EventList eventsData;  
  theAppBase.readEvents(eventsData, dataFileNames, withEvtWeight, 0, spinDensityHist::MAX_EVENTS);  
  
  EventList mcData; 
  theAppBase.readEvents(mcData, mcFileNames, false, 0, spinDensityHist::MAX_EVENTS);  

  std::shared_ptr<EvtDataBaseList> pbarpEventListPtr(new EvtDataBaseList(pbarpEnv::instance()));
  pbarpEventListPtr->read(eventsData, mcData);

  theLhPtr->setDataVec(pbarpEventListPtr->getDataVecs());
  theLhPtr->setMcVec(pbarpEventListPtr->getMcVecs());

  PwaFcnBase theFcn(theLhPtr, theFitParamBase, outputFileNameSuffix);

  std::shared_ptr<spinDensityHist> theSpinDensityHist(new spinDensityHist(theLhPtr, theStartparams));

  std::string serializationFileName = pbarpEnv::instance()->serializationFileName();
  std::ifstream serializationStream(serializationFileName.c_str());

  if(!serializationStream.is_open()){
     Warning << "Could not open serialization file." << endmsg;
  }
  else{
     boost::archive::text_iarchive boostInputArchive(serializationStream);
     std::shared_ptr<PwaCovMatrix> thePwaCovMatrix(new PwaCovMatrix);
     boostInputArchive >> *thePwaCovMatrix;
     theSpinDensityHist->SetCovarianceMatrix(thePwaCovMatrix);
  }
  theSpinDensityHist->Calculate();

  return 1;
 }


  EventList eventsData;  
  theAppBase.readEvents(eventsData, dataFileNames, withEvtWeight);  

  int ratioMcToData=theAppParams->ratioMcToData();
  int maxMcEvts=eventsData.size()*ratioMcToData;  
  EventList mcData; 
  theAppBase.readEvents(mcData, mcFileNames, false, 0, maxMcEvts-1);  


 if(mode == "server"){
   theAppBase.fixParams(upar,fixedParams); 

   std::shared_ptr<NetworkServer> theServer(new NetworkServer(theAppParams->serverPort(),
 								 theAppParams->noOfClients(),
 								 eventsData.size(),
 								 mcData.size()));

    PwaFcnServer theFcnServer(theLhPtr, theFitParamBase, theServer, outputFileNameSuffix);
    theServer->WaitForFirstClientLogin();

    FunctionMinimum min=theAppBase.migradDefault(theFcnServer, upar); 
 
    theServer->BroadcastClosingMessage();
    Info << "Closing server." << endmsg;

    theAppBase.printFitResult(min, theStartparams, std::cout, outputFileNameSuffix);

    return 1;
 }


 if(mode == "evoserver"){
   theAppBase.fixParams(upar,fixedParams); 

   std::shared_ptr<NetworkServer> theServer(new NetworkServer(theAppParams->serverPort(),
							      theAppParams->noOfClients(),
							      eventsData.size(),
							      mcData.size()));
   
   PwaFcnServer theFcnServer(theLhPtr, theFitParamBase, theServer, outputFileNameSuffix);
   theServer->WaitForFirstClientLogin();

   EvoMinimizer theEvoMinimizer(theFcnServer, upar, pbarpEnv::instance()->parser()->evoPopulation(),
				pbarpEnv::instance()->parser()->evoIterations());
   Info <<"start evolutionary minimizer "<< endmsg;
   std::vector<double> finalParamVec = theEvoMinimizer.Minimize();
   
   theServer->BroadcastClosingMessage();
   Info << "Closing server." << endmsg;
   
   fitParams finalFitParams=theStartparams;
   theFitParamBase->getFitParamVal(finalParamVec, finalFitParams);
   
   fitParams finalFitErrs=theErrorparams;
   
   std::ostringstream finalResultname;
   finalResultname << "finalResult" << outputFileNameSuffix << ".dat";

   std::ofstream theStream ( finalResultname.str().c_str() );
   theFitParamBase->dumpParams(theStream, finalFitParams, finalFitErrs);
   
   return 1;
 }


 std::shared_ptr<EvtDataBaseList> pbarpEventListPtr(new EvtDataBaseList(pbarpEnv::instance()));
 pbarpEventListPtr->read(eventsData, mcData);
 
 theLhPtr->setDataVec(pbarpEventListPtr->getDataVecs());
 theLhPtr->setMcVec(pbarpEventListPtr->getMcVecs());
 
 PwaFcnBase theFcn(theLhPtr, theFitParamBase, outputFileNameSuffix); 
 Info << "\nThe parameter values are: " << "\n" << endmsg;
 theFitParamBase->printParams(theStartparams);
 
 Info << "\nThe parameter errors are: " << "\n" << endmsg;
 theFitParamBase->printParams(theErrorparams);

 
 if (mode=="qaMode"){
   double evtWeightSumData = pbarpEventListPtr->NoOfWeightedDataEvts();
   theAppBase.qaMode(theStartparams, evtWeightSumData, noOfFreeFitParams );
   pbarpHist theHist(theLhPtr, theStartparams);     
     end= clock();
   double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
   Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;
    
    return 1;    
  }
  

  if (mode=="pwa"){

    bool cacheAmps = theAppParams->cacheAmps();
    Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
    if (cacheAmps) theLhPtr->cacheAmplitudes();

    theAppBase.fixParams(upar,fixedParams); 
     FunctionMinimum min=theAppBase.migradDefault(theFcn, upar);
    theAppBase.printFitResult(min, theStartparams, std::cout, outputFileNameSuffix); 

    return 1;
 }

  if (mode=="evo"){

    bool cacheAmps = theAppParams->cacheAmps();
    Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
    if (cacheAmps) theLhPtr->cacheAmplitudes();

    theAppBase.fixParams(upar,fixedParams); 

    EvoMinimizer theEvoMinimizer(theFcn, upar, pbarpEnv::instance()->parser()->evoPopulation(),
				 pbarpEnv::instance()->parser()->evoIterations());
    Info <<"start evolutionary minimizer "<< endmsg;
    std::vector<double> finalParamVec = theEvoMinimizer.Minimize();

    fitParams finalFitParams=theStartparams;
    theFitParamBase->getFitParamVal(finalParamVec, finalFitParams);

    fitParams finalFitErrs=theErrorparams;

    std::ostringstream finalResultname;
    finalResultname << "finalResult" << outputFileNameSuffix << ".dat";

    std::ofstream theStream ( finalResultname.str().c_str() );
    theFitParamBase->dumpParams(theStream, finalFitParams, finalFitErrs);

    return 1;
  }

 return 1;
}     
 
  
