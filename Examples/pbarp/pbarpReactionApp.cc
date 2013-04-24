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
#include <boost/shared_ptr.hpp>

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
#include "PwaUtils/PwaGen.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/ErrLogUtils.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "pbarpUtils/pbarpReaction.hh"
#include "pbarpUtils/pbarpBaseLh.hh"
#include "pbarpUtils/pbarpHeliLh.hh"
#include "pbarpUtils/pbarpCanoLh.hh"

#include "Event/EventReaderDefault.hh"

#include "PwaUtils/EvtDataBaseList.hh"
#include "pbarpUtils/pbarpHist.hh"
#include "pbarpUtils/spinDensityHist.hh"
#include "Event/Event.hh"
#include "Event/EventList.hh"
#include "PwaUtils/NetworkServer.hh"
#include "PwaUtils/NetworkClient.hh"
#include "PwaUtils/WelcomeScreen.hh"

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

  boost::shared_ptr<pbarpReaction> thepbarpReaction=pbarpEnv::instance()->reaction();

  thepbarpReaction->print(std::cout);


  std::string mode=theAppParams->mode();

  boost::shared_ptr<FitParamsBase> theFitParamBase=boost::shared_ptr<FitParamsBase>(new FitParamsBase());


  std::string prodFormalism=theAppParams->productionFormalism();
  boost::shared_ptr<AbsLh> theLhPtr;
  if(prodFormalism=="Cano") theLhPtr=boost::shared_ptr<AbsLh>(new pbarpCanoLh());
  else if(prodFormalism=="Heli") theLhPtr=boost::shared_ptr<AbsLh>(new pbarpHeliLh());
  else {
    Alert << "prodFormalism\t" << prodFormalism << "\tdoesn't exist!!!" << endmsg;
    exit(1);
  }

  if (mode=="dumpDefaultParams"){
    fitParams defaultVal;
    fitParams defaultErr;
    theLhPtr->getDefaultParams(defaultVal, defaultErr);

    std::stringstream defaultparamsname;
    defaultparamsname << "defaultparams" << pbarpEnv::instance()->outputFileNameSuffix() << ".dat";
    std::ofstream theStreamDefault ( defaultparamsname.str().c_str() );
    
    theFitParamBase->dumpParams(theStreamDefault, defaultVal, defaultErr);
    return 1;
  }


  std::string paramStreamerPath=theAppParams->fitParamFile();
  std::string outputFileNameSuffix= pbarpEnv::instance()->outputFileNameSuffix();
  StreamFitParmsBase theParamStreamer(paramStreamerPath, theLhPtr);
  fitParams theStartparams=theParamStreamer.getFitParamVal();
  fitParams theErrorparams=theParamStreamer.getFitParamErr();

  if (mode=="gen"){
    boost::shared_ptr<PwaGen> pwaGenPtr(new PwaGen(pbarpEnv::instance()));
    pwaGenPtr->generate(theLhPtr, theStartparams);
    theFitParamBase->printParams(theStartparams);
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

  
  int noFinalStateParticles=pbarpEnv::instance()->noFinalStateParticles();  



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
  
  
  EventReaderDefault eventReaderDataClient(dataFileNames, noFinalStateParticles, 0, withEvtWeight);
  eventReaderDataClient.setUnit(theAppParams->unitInFile());
  eventReaderDataClient.setOrder(theAppParams->orderInFile());
  
  EventList* eventsDataClient=new EventList();
  eventReaderDataClient.fill(*eventsDataClient, theClient.GetEventLimits()[0], theClient.GetEventLimits()[1]);
  
  eventsDataClient->rewind();  
  Info  << "\nFile has " << eventsDataClient->size() << " events. Each event has "
        <<  eventsDataClient->nextEvent()->size() << " final state particles.\n" ;  // << endmsg;
  eventsDataClient->rewind();
  
  
  
  EventReaderDefault eventReaderMcClient(mcFileNames, noFinalStateParticles, 0, false);
  eventReaderMcClient.setUnit(theAppParams->unitInFile());
  eventReaderMcClient.setOrder(theAppParams->orderInFile());
  
  
  EventList* mcDataClient=new EventList();
  eventReaderMcClient.fill(*mcDataClient, theClient.GetEventLimits()[2], theClient.GetEventLimits()[3]);
  Info  << "\nFile has " << mcDataClient->size() << " events. Each event has "
        <<  mcDataClient->nextEvent()->size() << " final state particles.\n" ;  // << endmsg;
  mcDataClient->rewind();
  
  boost::shared_ptr<EvtDataBaseList> pbarpEventListPtr(new EvtDataBaseList(pbarpEnv::instance()));
  pbarpEventListPtr->read(*eventsDataClient, *mcDataClient);

  delete eventsDataClient;
  delete mcDataClient;
  
  theLhPtr->setDataVec(pbarpEventListPtr->getDataVecs());
  theLhPtr->setMcVec(pbarpEventListPtr->getMcVecs());
  
  PwaFcnBase theFcn(theLhPtr, theFitParamBase, outputFileNameSuffix);
  
  while(true){
    
    if(!theClient.WaitForParams())
      return 0;
    
    const std::vector<double> currentParamVec=theClient.GetParams();
    fitParams currentFitParams=theStartparams;
    theFitParamBase->getFitParamVal(currentParamVec, currentFitParams);
    
    LHData theLHData;
    theLhPtr->calcLogLhDataClient(currentFitParams, theLHData);
    
    if(!theClient.SendLH(theLHData.logLH_data, theLHData.weightSum, theLHData.LH_mc))
      return 0;
  }
  return 1;
 }


  EventReaderDefault eventReaderData(dataFileNames, noFinalStateParticles, 0, withEvtWeight);
  eventReaderData.setUnit(theAppParams->unitInFile());
  eventReaderData.setOrder(theAppParams->orderInFile());

  EventList eventsData;
  eventReaderData.fill(eventsData);
  
  Info  << "\nFile has " << eventsData.size() << " events. Each event has "
	<<  eventsData.nextEvent()->size() << " final state particles.\n" ;  // << endmsg;
  eventsData.rewind();

  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = eventsData.nextEvent()) != 0 && evtCount < 10) {
    Info        << "\n";
    for(int i=0; i<noFinalStateParticles; ++i){
      Info        << (*anEvent->p4(i)) << "\tm = " << anEvent->p4(i)->Mass() << "\n";
    }
    Info        << "\n" << endmsg;
    ;  // << endmsg;
    ++evtCount;
  }
  eventsData.rewind();

  EventReaderDefault eventReaderMc(mcFileNames, noFinalStateParticles, 0, false);
  eventReaderMc.setUnit(theAppParams->unitInFile());
  eventReaderMc.setOrder(theAppParams->orderInFile());

  int ratioMcToData=theAppParams->ratioMcToData();
  int maxMcEvts=eventsData.size()*ratioMcToData;

  EventList mcData;
  eventReaderMc.fill(mcData, 0, maxMcEvts-1);
  Info  << "\nFile has " << mcData.size() << " events. Each event has "
        <<  mcData.nextEvent()->size() << " final state particles.\n" ;  // << endmsg;
  mcData.rewind();

  evtCount = 0;
  while ((anEvent = mcData.nextEvent()) != 0 && evtCount < 10) {
    Info        << "\n";
    for(int i=0; i<noFinalStateParticles; ++i){
      Info        << (*anEvent->p4(i)) << "\tm = " << anEvent->p4(i)->Mass() << "\n";
    }
    Info        << "\n" << endmsg;
    ;  // << endmsg;
    ++evtCount;
  }
  mcData.rewind();

 if(mode == "server"){

    std::vector<std::string>::const_iterator itFix;
    for (itFix=fixedParams.begin(); itFix!=fixedParams.end(); ++itFix){
       upar.Fix( (*itFix) );
    }

    boost::shared_ptr<NetworkServer> theServer(new NetworkServer(theAppParams->serverPort(),
 								 theAppParams->noOfClients(),
 								 eventsData.size(),
 								 mcData.size()));

    PwaFcnServer theFcnServer(theLhPtr, theFitParamBase, theServer, outputFileNameSuffix);
    theServer->WaitForFirstClientLogin();

    MnMigrad migrad(theFcnServer, upar);
    Info <<"start migrad "<< endmsg;
    FunctionMinimum min = migrad();

    if(!min.IsValid()) {
      //try with higher strategy
      Info <<"FM is invalid, try with strategy = 2."<< endmsg;
      MnMigrad migrad2(theFcnServer, min.UserState(), MnStrategy(2));
      min = migrad2();
    }

    theServer->SendClosingMessage();
    Info << "Closing server." << endmsg;


    std::cout << "\n\n********************** Final fit parameters *************************" << std::endl;
    std::cout << "\n" << min.UserParameters() << std::endl;
    std::cout << "\n\n**************** Minuit FunctionMinimum information ******************" << std::endl;
    if(min.IsValid())             std::cout << "\n Function minimum is valid." << std::endl;
    else                          std::cout << "\n WARNING: Function minimum is invalid!" << std::endl;
    if(min.HasValidCovariance())  std::cout << "\n Covariance matrix is valid." << std::endl;
    else                          std::cout << "\n WARNING: Covariance matrix is invalid!" << std::endl;
    std::cout <<"\n Final LH: "<< std::setprecision(10) << min.Fval() << "\n" << std::endl;
    std::cout <<" # of function calls: " << min.NFcn() << std::endl;
    std::cout <<" minimum edm: " << std::setprecision(10) << min.Edm()<<std::endl;
    if(!min.HasValidParameters()) std::cout << " hasValidParameters() returned FALSE" << std::endl;
    if(!min.HasAccurateCovar())   std::cout << " hasAccurateCovar() returned FALSE" << std::endl;
    if(!min.HasPosDefCovar())     std::cout << " hasPosDefCovar() returned FALSE" << std::endl;
    if(!min.HasMadePosDefCovar()) std::cout << " hasMadePosDefCovar() returned FALSE" << std::endl;
    if(!min.HasCovariance())      std::cout << " hasCovariance() returned FALSE" << std::endl;
    if(min.HasReachedCallLimit()) std::cout << " hasReachedCallLimit() returned TRUE" << std::endl;
    if(min.IsAboveMaxEdm())       std::cout << " isAboveMaxEdm() returned TRUE" << std::endl;
    if(min.HesseFailed())         std::cout << " hesseFailed() returned TRUE" << std::endl;
    std::cout << std::endl;



    MnUserParameters finalUsrParameters=min.UserParameters();
    const std::vector<double> finalParamVec=finalUsrParameters.Params();
    fitParams finalFitParams=theStartparams;
    theFitParamBase->getFitParamVal(finalParamVec, finalFitParams);

    const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
    fitParams finalFitErrs=theErrorparams;
    theFitParamBase->getFitParamVal(finalParamErrorVec, finalFitErrs);

    std::ostringstream finalResultname;
    finalResultname << "finalResult" << outputFileNameSuffix << ".dat";

    std::ofstream theStream ( finalResultname.str().c_str() );
    theFitParamBase->dumpParams(theStream, finalFitParams, finalFitErrs);

    MnUserCovariance theCovMatrix = min.UserCovariance();
    std::ostringstream serializationFileName;
    serializationFileName << "serializedOutput" << outputFileNameSuffix << ".dat";
    std::ofstream serializationStream(serializationFileName.str().c_str());
    boost::archive::text_oarchive boostOutputArchive(serializationStream);

    if(min.IsValid()){
       PwaCovMatrix thePwaCovMatrix(theCovMatrix, finalUsrParameters, finalFitParams);
       boostOutputArchive << thePwaCovMatrix;
    }
    return 1;
 }


 boost::shared_ptr<EvtDataBaseList> pbarpEventListPtr(new EvtDataBaseList(pbarpEnv::instance()));
 pbarpEventListPtr->read(eventsData, mcData);
 
 theLhPtr->setDataVec(pbarpEventListPtr->getDataVecs());
 theLhPtr->setMcVec(pbarpEventListPtr->getMcVecs());
 
 PwaFcnBase theFcn(theLhPtr, theFitParamBase, outputFileNameSuffix); 
 Info << "\nThe parameter values are: " << "\n" << endmsg;
 theFitParamBase->printParams(theStartparams);
 
 Info << "\nThe parameter errors are: " << "\n" << endmsg;
 theFitParamBase->printParams(theErrorparams);

 
 if (mode=="qaMode"){

    double theLh=theLhPtr->calcLogLh(theStartparams);
    Info <<"theLh = "<< theLh << endmsg;

    pbarpHist theHist(theLhPtr, theStartparams);

    double evtWeightSumData = pbarpEventListPtr->NoOfWeightedDataEvts();
    double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
    double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
    double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/(evtWeightSumData-noOfFreeFitParams-1);
    
    boost::shared_ptr<WaveContribution> theWaveContribution;
    if(pbarpEnv::instance()->parser()->calcContributionError()){
       std::string serializationFileName = pbarpEnv::instance()->serializationFileName();
       std::ifstream serializationStream(serializationFileName.c_str());

       if(!serializationStream.is_open()){
	  Alert << "Could not open serialization file." << endmsg;
	  return 0;
       }

       boost::archive::text_iarchive boostInputArchive(serializationStream);

       boost::shared_ptr<PwaCovMatrix> thePwaCovMatrix(new PwaCovMatrix);
       boostInputArchive >> *thePwaCovMatrix;
       theWaveContribution = boost::shared_ptr<WaveContribution>
	  (new WaveContribution(theLhPtr, theStartparams, thePwaCovMatrix));
    }
    else{
       theWaveContribution = boost::shared_ptr<WaveContribution>
	  (new WaveContribution(theLhPtr, theStartparams));
    }

    std::pair<double, double> contValue = theWaveContribution->CalcContribution();

    Info << "noOfFreeFitParams:\t" <<noOfFreeFitParams;
    Info << "evtWeightSumData:\t" <<evtWeightSumData; 
    Info << "BIC:\t" << BICcriterion << endmsg;
    Info << "AIC:\t" << AICcriterion << endmsg;
    Info << "AICc:\t" << AICccriterion << endmsg;
    Info << "Selected wave contribution:\t" << contValue.first
	 << " +- " << contValue.second << endmsg;

    std::ostringstream qaSummaryFileName;
    qaSummaryFileName << "qaSummary" << outputFileNameSuffix << ".dat";

    std::ofstream theQaStream ( qaSummaryFileName.str().c_str() );
    theQaStream << "BIC\t" << BICcriterion << "\n";
    theQaStream << "AICa\t" << AICcriterion << "\n";
    theQaStream << "AICc\t" << AICccriterion << "\n";
    theQaStream << "logLh\t" << theLh << "\n";
    theQaStream << "free parameter\t" << noOfFreeFitParams << "\n";
    theQaStream << "Selected wave contribution\t" << contValue.first
		<< " +- " << contValue.second <<  "\n";
    theQaStream.close();
    
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
    

    std::cout << "\n\n********************** Final fit parameters *************************" << std::endl;
    std::cout << "\n" << min.UserParameters() << std::endl;
    std::cout << "\n\n**************** Minuit FunctionMinimum information ******************" << std::endl;
    if(min.IsValid())             std::cout << "\n Function minimum is valid." << std::endl;
    else                          std::cout << "\n WARNING: Function minimum is invalid!" << std::endl;
    if(min.HasValidCovariance())  std::cout << "\n Covariance matrix is valid." << std::endl;
    else                          std::cout << "\n WARNING: Covariance matrix is invalid!" << std::endl;
    std::cout <<"\n Final LH: "<< std::setprecision(10) << min.Fval() << "\n" << std::endl;
    std::cout <<" # of function calls: " << min.NFcn() << std::endl;
    std::cout <<" minimum edm: " << std::setprecision(10) << min.Edm()<<std::endl;    
    if(!min.HasValidParameters()) std::cout << " hasValidParameters() returned FALSE" << std::endl;
    if(!min.HasAccurateCovar())   std::cout << " hasAccurateCovar() returned FALSE" << std::endl;
    if(!min.HasPosDefCovar())     std::cout << " hasPosDefCovar() returned FALSE" << std::endl;
    if(!min.HasMadePosDefCovar()) std::cout << " hasMadePosDefCovar() returned FALSE" << std::endl;
    if(!min.HasCovariance())      std::cout << " hasCovariance() returned FALSE" << std::endl;
    if(min.HasReachedCallLimit()) std::cout << " hasReachedCallLimit() returned TRUE" << std::endl;
    if(min.IsAboveMaxEdm())       std::cout << " isAboveMaxEdm() returned TRUE" << std::endl;
    if(min.HesseFailed())         std::cout << " hesseFailed() returned TRUE" << std::endl;
    std::cout << std::endl;


    MnUserParameters finalUsrParameters=min.UserParameters();
    const std::vector<double> finalParamVec=finalUsrParameters.Params();
    fitParams finalFitParams=theStartparams;
    theFitParamBase->getFitParamVal(finalParamVec, finalFitParams);

    const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
    fitParams finalFitErrs=theErrorparams;
    theFitParamBase->getFitParamVal(finalParamErrorVec, finalFitErrs);
    
    std::ostringstream finalResultname;
    finalResultname << "finalResult" << outputFileNameSuffix << ".dat";

    std::ofstream theStream ( finalResultname.str().c_str() );
    theFitParamBase->dumpParams(theStream, finalFitParams, finalFitErrs);
    
    MnUserCovariance theCovMatrix = min.UserCovariance();

    std::ostringstream serializationFileName;
    serializationFileName << "serializedOutput" << outputFileNameSuffix << ".dat";
    std::ofstream serializationStream(serializationFileName.str().c_str());
    boost::archive::text_oarchive boostOutputArchive(serializationStream);

    if(min.HasValidCovariance()){
       PwaCovMatrix thePwaCovMatrix(theCovMatrix, finalUsrParameters, finalFitParams);
       boostOutputArchive << thePwaCovMatrix;
    }

    return 1;
 }

 if(mode == "spinDensity"){
 
    std::string serializationFileName = pbarpEnv::instance()->serializationFileName();
    std::ifstream serializationStream(serializationFileName.c_str());

    if(!serializationStream.is_open()){
       Alert << "Could not open serialization file." << endmsg;
       return 0;
    }

    boost::archive::text_iarchive boostInputArchive(serializationStream);

    PwaCovMatrix thePwaCovMatrix;
    boostInputArchive >> thePwaCovMatrix;

    spinDensityHist theSpinDensityHists(theLhPtr, theStartparams, thePwaCovMatrix);
    return 1;
 }

 return 1;
}     
 
  
