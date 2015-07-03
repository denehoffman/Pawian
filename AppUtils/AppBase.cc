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

// AppBase class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <string>
#include <memory>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AppUtils/AppBase.hh"

#include "AppUtils/LhFactory.hh"

#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtWeightList.hh"
#include "PwaUtils/PwaGen.hh"
#include "PwaUtils/AbsHist.hh"
#include "PwaUtils/WaveContribution.hh"
#include "PwaUtils/NetworkServer.hh"
#include "PwaUtils/NetworkClient.hh"
#include "PwaUtils/EvtDataBaseList.hh"

#include "FitParams/AbsPawianParamStreamer.hh"
#include "FitParams/PwaCovMatrix.hh"

#include "ConfigParser/ParserBase.hh"
#include "ConfigParser/pbarpParser.hh"
#include "ConfigParser/epemParser.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Event/Event.hh"
#include "Event/EventReaderDefault.hh"

#include "MinFunctions/PwaFcnBase.hh"
#include "MinFunctions/PwaFcnServer.hh"
#include "MinFunctions/AbsPawianMinimizer.hh"
#include "MinFunctions/EvoMinimizer.hh"
#include "MinFunctions/MinuitMinimizer.hh"

#include "FitParams/AbsPawianParameters.hh"
#include "FitParams/ParamDepHandler.hh"

#include "pbarpUtils/PbarpChannelEnv.hh"
#include "epemUtils/EpemChannelEnv.hh"

AppBase::AppBase()
{
   // Check mode for coupled channel analyses
   std::string mode = GlobalEnv::instance()->parser()->mode();
   if((GlobalEnv::instance()->NoChannels() > 1) &&
         ((mode == "spinDensity") || (mode == "pwa") || (mode == "qaMode"))){
      Alert << "Selected mode does not support coupled channel analyses!" << endmsg;
   }
}

AppBase::~AppBase()
{
}

void AppBase::dumpDefaultParams(){


    std::stringstream defaultparamsname;
    defaultparamsname << "defaultparams" << GlobalEnv::instance()->outputFileNameSuffix() << ".dat";
    std::ofstream theStreamDefault ( defaultparamsname.str().c_str() );

    std::shared_ptr<AbsPawianParameters> defaultParams=GlobalEnv::instance()->defaultPawianParams();
    defaultParams->print(theStreamDefault);
}

void AppBase::generate(std::shared_ptr<AbsPawianParameters> theParams){
    std::shared_ptr<PwaGen> pwaGenPtr(new PwaGen());
    theParams->print(std::cout);
    pwaGenPtr->generate(GlobalEnv::instance()->Channel()->Lh(), theParams);
}

void AppBase::readEvents(EventList& theEventList, std::vector<std::string>& fileNames, ChannelID channelID, bool withEvtWeight, int evtStart, int evtStop){
  int noFinalStateParticles=GlobalEnv::instance()->Channel(channelID)->noFinalStateParticles();
  std::vector< std::shared_ptr<MassRangeCut> > massRangeCuts=GlobalEnv::instance()->Channel(channelID)->massRangeCuts();
  EventReaderDefault eventReader(fileNames, noFinalStateParticles, 0, withEvtWeight);
  eventReader.setUnit(GlobalEnv::instance()->Channel(channelID)->parser()->unitInFile());
  eventReader.setOrder(GlobalEnv::instance()->Channel(channelID)->parser()->orderInFile());

  if(GlobalEnv::instance()->Channel(channelID)->useMassRange()){
   eventReader.setMassRange(massRangeCuts);
  }

  eventReader.fill(theEventList, evtStart, evtStop);

  Info  << "\nFile has " << theEventList.size() << " events. Each event has "
        <<  theEventList.nextEvent()->size() << " final state particles.\n" ;  // << endmsg;
  theEventList.rewind();

  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = theEventList.nextEvent()) != 0 && evtCount < 10) {
    Info        << "\n";
    for(int i=0; i<noFinalStateParticles; ++i){
      Info        << (*anEvent->p4(i)) << "\tm = " << anEvent->p4(i)->Mass() << "\n";
    }
    Info        << "\n" << endmsg;
    ++evtCount;
  }
  theEventList.rewind();
}

void AppBase::createLhObjects(){

  ChannelEnvList channelEnvs=GlobalEnv::instance()->ChannelEnvs();
  for(auto it=channelEnvs.begin();it!=channelEnvs.end();++it){
    std::string prodFormalism=(*it).first->parser()->productionFormalism();
    std::shared_ptr<AbsLh> theLhPtr=LhFactory::instance()->getLh((*it).first->channelType(),(*it).first->channelID(), prodFormalism );
    (*it).first->SetLh(theLhPtr);
  }
  
  // Generate the full parameter set using the likelihood list
  //  GlobalEnv::instance()->CreateDefaultParameterSet();
}

void AppBase::qaMode(std::shared_ptr<AbsPawianParameters> startParams, double evtWeightSumData){
  int noOfFreeFitParams=startParams->VariableParameters();
  if(evtWeightSumData<=(noOfFreeFitParams+1)){
    Warning << "number of data events less or equal to the number of free parameters!!!"
	    << "\n be careful with the fit result!!!"
	    << endmsg; 
      }  
  
  double theLh=GlobalEnv::instance()->Channel()->Lh()->calcLogLh(startParams);
  double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
  double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
  double AICccDenom = evtWeightSumData-noOfFreeFitParams-1;
  if ( fabs(AICccDenom) < 1.e-10) AICccDenom=1.e-10;
  double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/AICccDenom;


  std::shared_ptr<WaveContribution> theWaveContribution;
  if(GlobalEnv::instance()->parser()->calcContributionError()){
    std::string serializationFileName = GlobalEnv::instance()->serializationFileName();
    std::ifstream serializationStream(serializationFileName.c_str());
    
    if(!serializationStream.is_open()){
      Alert << "Could not open serialization file." << endmsg;
      exit(0);
    }
    boost::archive::text_iarchive boostInputArchive(serializationStream);
    
    std::shared_ptr<PwaCovMatrix> thePwaCovMatrix(new PwaCovMatrix);
    boostInputArchive >> *thePwaCovMatrix;
    theWaveContribution = std::shared_ptr<WaveContribution>
      (new WaveContribution(GlobalEnv::instance()->Channel()->Lh(), startParams, thePwaCovMatrix));
  }
  else{
    theWaveContribution = std::shared_ptr<WaveContribution>
      (new WaveContribution(GlobalEnv::instance()->Channel()->Lh(), startParams));
  }
  std::pair<double, double> contValue = theWaveContribution->CalcContribution();
  std::vector<std::pair<std::string,std::pair<double,double>>> singleContValues = theWaveContribution->CalcSingleContributions();

  std::ostringstream qaSummaryFileName;
  std::string outputFileNameSuffix= GlobalEnv::instance()->outputFileNameSuffix();
  qaSummaryFileName << "qaSummary" << outputFileNameSuffix << ".dat";
  std::ofstream theQaStream ( qaSummaryFileName.str().c_str() );

  Info        << "logLh\t" << theLh;
  theQaStream << "logLh\t" << theLh << "\n";

  Info        << "noOfFreeFitParams:\t" << noOfFreeFitParams;
  theQaStream << "noOfFreeFitParams:\t" << noOfFreeFitParams << "\n";

  Info        << "BIC:\t" << BICcriterion;
  theQaStream << "BIC:\t" << BICcriterion << "\n";

  Info        << "AICa:\t" << AICcriterion;
  theQaStream << "AICa:\t" << AICcriterion << "\n";

  Info        << "AICc:\t" << AICccriterion;
  theQaStream << "AICc:\t" << AICccriterion << "\n";

  Info        << "No of data events without weight " << GlobalEnv::instance()->Channel()->Lh()->getDataVec().size();
  theQaStream << "No of data events without weight " << GlobalEnv::instance()->Channel()->Lh()->getDataVec().size() << "\n";

  Info        << "No of data events with weight " << evtWeightSumData;
  theQaStream << "No of data events with weight " << evtWeightSumData << "\n";

  Info        << "No of MC events " << GlobalEnv::instance()->Channel()->Lh()->getMcVec().size();
  theQaStream << "No of MC events " << GlobalEnv::instance()->Channel()->Lh()->getMcVec().size() << "\n";

  double scaleFactor = evtWeightSumData/GlobalEnv::instance()->Channel()->Lh()->getMcVec().size();
  Info        << "scaling factor " << scaleFactor;
  theQaStream << "scaling factor " << scaleFactor << "\n";

  Info        << "no of fitted events with scaling factor: " << contValue.first*scaleFactor;
  theQaStream << "no of fitted events with scaling factor: " << contValue.first*scaleFactor << "\n";

  Info        << "Selected wave contribution:\t" << contValue.first << " +- " << contValue.second;
  theQaStream << "Selected wave contribution:\t" << contValue.first << " +- " << contValue.second <<  "\n";

  std::vector<std::pair<std::string,std::pair<double,double>>>::iterator it;
  for(it=singleContValues.begin(); it!=singleContValues.end(); ++it) {
    Info        << "Single wave contribution " << (*it).first << "\t" << (*it).second.first << " +- " << (*it).second.second;
    theQaStream << "Single wave contribution " << (*it).first << "\t" << (*it).second.first << " +- " << (*it).second.second <<  "\n";
  }

  theQaStream.close();

  std::shared_ptr<AbsHist> histPtr = GlobalEnv::instance()->Channel()->CreateHistInstance();
  histPtr->fillFromLhData(GlobalEnv::instance()->Channel()->Lh(), startParams);
}

void AppBase::qaModeSimple(EventList& dataEventList, EventList& mcEventList, std::shared_ptr<AbsPawianParameters> startParams){
  std::shared_ptr<EvtDataBaseList> evtDataBaseList(new EvtDataBaseList(0));
  int noOfFreeFitParams=startParams->VariableParameters();
  std::shared_ptr<AbsLh> absLh=GlobalEnv::instance()->Channel()->Lh();
  LHData theLHData;
  std::shared_ptr<WaveContribution> theWaveContribution(new WaveContribution(GlobalEnv::instance()->Channel()->Lh(), startParams));

  std::shared_ptr<AbsPawianParameters> currentParams = std::shared_ptr<AbsPawianParameters>(startParams->Clone());

  std::ostringstream qaSummaryFileName;
  std::string outputFileNameSuffix= GlobalEnv::instance()->outputFileNameSuffix();
  qaSummaryFileName << "qaSummarySimple" << outputFileNameSuffix << ".dat";
  std::ofstream theQaStream ( qaSummaryFileName.str().c_str() );

  for(int i=-1; i<static_cast<int>(theWaveContribution->NoOfContributions());i++){

    std::string contributionName="";

    if(i!=-1){
      contributionName = theWaveContribution->GetContributionName(i);
      currentParams = theWaveContribution->GetParametersForContribution(i);
    }
 
    std::shared_ptr<AbsHist> histPtr = GlobalEnv::instance()->Channel()->CreateHistInstance(contributionName);

    absLh->updateFitParams(currentParams);
  
    //loop over data events
    Event* anEvent;
    int evtCount = 0;

    double integralDataWoWeight=(double) dataEventList.size();
    double evtWeightSumData=0.;
    
    dataEventList.rewind();
    while ((anEvent = dataEventList.nextEvent())){
      EvtData* currentDataEvt=evtDataBaseList->convertEvent(anEvent, evtCount);
      absLh->addDataToLogLh(currentDataEvt, currentParams, theLHData);
      histPtr->fillEvt(currentDataEvt, currentDataEvt->evtWeight, "data");
      evtWeightSumData += currentDataEvt->evtWeight;
      delete currentDataEvt;
      evtCount++;
      if (evtCount%1000 == 0) Info << evtCount << " data events calculated" << endmsg;
    }

    //loop over mc events
    int evtCountMc = 0;

    double integralFitWeight=0.;

    mcEventList.rewind();
    while ((anEvent = mcEventList.nextEvent())){
      EvtData* currentMcEvt=evtDataBaseList->convertEvent(anEvent, evtCount);
      double currentIntensity=absLh->addMcToLogLh(currentMcEvt,currentParams, theLHData);
      histPtr->fillEvt(currentMcEvt, 1., "mc");
      histPtr->fillEvt(currentMcEvt, currentIntensity, "fit");

      integralFitWeight+=currentIntensity;

      delete currentMcEvt;
      evtCount++;
      evtCountMc++;
      if (evtCountMc%1000 == 0) Info << evtCountMc << " MC events calculated" << endmsg ;
    }

    double scaleFactor=theLHData.weightSum/theLHData.num_mc;
    histPtr->scaleFitHists(scaleFactor);

    // if(i!=-1)
    //   continue;
    if(i==-1){
      double theLh=absLh->mergeLogLhData(theLHData);
      //    double evtWeightSumData=theLHData.weightSum;
      double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
      double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
      double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/(evtWeightSumData-noOfFreeFitParams-1);
      //    double integralDataWoWeight=(double) dataEventList.size();
      
      // std::ostringstream qaSummaryFileName;
      // std::string outputFileNameSuffix= GlobalEnv::instance()->outputFileNameSuffix();
      // qaSummaryFileName << "qaSummarySimple" << outputFileNameSuffix << ".dat";
      // std::ofstream theQaStream ( qaSummaryFileName.str().c_str() );
      
      Info        << "logLh\t" << theLh;
      theQaStream << "logLh\t" << theLh << "\n";
      
      Info        << "noOfFreeFitParams:\t" << noOfFreeFitParams;
      theQaStream << "noOfFreeFitParams\t" << noOfFreeFitParams << "\n";
      
      Info        << "BIC:\t" << BICcriterion;
      theQaStream << "BIC:\t" << BICcriterion << "\n";
      
      Info        << "AICa:\t" << AICcriterion;
      theQaStream << "AICa:\t" << AICcriterion << "\n";
      
      Info        << "AICc:\t" << AICccriterion;
      theQaStream << "AICc:\t" << AICccriterion << "\n";
      
      Info        << "No of data events without weight " << integralDataWoWeight;
      theQaStream << "No of data events without weight " << integralDataWoWeight << "\n";
      
      Info        << "No of data events with weight " << evtWeightSumData;
      theQaStream << "No of data events with weight " << evtWeightSumData << "\n";
      
      Info        << "No of MC events " << theLHData.num_mc;
      theQaStream << "No of MC events " << theLHData.num_mc << "\n";
      
      Info        << "scaling factor " << scaleFactor;
      theQaStream << "scaling factor " << scaleFactor << "\n";
      
      Info        << "no of fitted events with scaling factor: " << integralFitWeight*scaleFactor;
      theQaStream << "no of fitted events with scaling factor: " << integralFitWeight*scaleFactor << "\n";
    }
    else{ //i>-1
      Info        << "contribution no " << i;
      Info        << "No of data events without weight " << integralDataWoWeight;
      Info        << "No of data events with weight " << evtWeightSumData;
      Info        << "No of MC events " << theLHData.num_mc;
      Info        << "no of fitted events with scaling factor: " << integralFitWeight*scaleFactor;
      theQaStream << "contribution no " << i << "\n";
      theQaStream << "no of fitted events with scaling factor: " << integralFitWeight*scaleFactor << "\n";
    }
  } // loop over contributions
  theQaStream.close();
}

void AppBase::plotMode(EventList& dataEventList, EventList& mcEventList, std::shared_ptr<EvtDataBaseList> evtDataBaseList){

  std::shared_ptr<AbsHist> histPtr = GlobalEnv::instance()->Channel()->CreateHistInstance();

  //loop over data events
  Event* anEvent;
  int evtCount = 0;
  double evtWeightSumData=0.;
  dataEventList.rewind();
  while ((anEvent = dataEventList.nextEvent())){
    // EvtData* currentDataEvt=evtDataBaseList->convertEvent(anEvent, evtCount);
    EvtData* currentDataEvt=evtDataBaseList->convertEvent(anEvent);
    histPtr->fillEvt(currentDataEvt, currentDataEvt->evtWeight, "data");
    evtWeightSumData+=currentDataEvt->evtWeight;
    delete currentDataEvt;
    evtCount++;
    if (evtCount%1000 == 0) Info << evtCount << " data events calculated" << endmsg;
  }

  //loop over mc events
  int evtCountMc = 0;
  mcEventList.rewind();
  while ((anEvent = mcEventList.nextEvent())){
    // EvtData* currentMcEvt=evtDataBaseList->convertEvent(anEvent, evtCount);
    EvtData* currentMcEvt=evtDataBaseList->convertEvent(anEvent);
    histPtr->fillEvt(currentMcEvt, 1., "fit");
    histPtr->fillEvt(currentMcEvt, 1., "mc");
    delete currentMcEvt;
    evtCount++;
    evtCountMc++;
    if (evtCountMc%1000 == 0) Info << evtCountMc << " MC events calculated" << endmsg ;
  }

  double histScaleFactor=evtWeightSumData/evtCountMc;
  histPtr->scaleFitHists(histScaleFactor);
}

std::shared_ptr<AbsPawianParameters> AppBase::streamPawianParams(){
  std::string paramStreamerPath=GlobalEnv::instance()->parser()->fitParamFile();
  AbsPawianParamStreamer thePawianStreamer(paramStreamerPath);
  return thePawianStreamer.paramList();
}

void AppBase::fixParams(std::shared_ptr<AbsPawianParameters> upar, std::vector<std::string> fixedParams, bool mustMatch){

  // Evaluate parameter dependencies and add fixes
  ParamDepHandler::instance()->Fill(GlobalEnv::instance()->parser()->parameterDependencies(), upar);
  std::vector<std::string> dependentParameters = ParamDepHandler::instance()->DependentParameterNames();
  // fixedParams.insert(fixedParams.end(), dependentParameters.begin(), dependentParameters.end());

  std::vector<std::string>::const_iterator strIt;
  for(strIt=dependentParameters.begin(); strIt!=dependentParameters.end(); ++strIt){
    std::vector<unsigned int> depIdVec=ParamDepHandler::instance()->RefIds(*strIt);
    upar->SetDependencies( *strIt, depIdVec);
  }

  // Always fix the primary channel's scaling parameters
  //  std::string fixedScaleParam = GlobalEnv::instance()->Channel()->Lh()->getChannelScaleParam() + "Other";
  std::string fixedScaleParam = GlobalEnv::instance()->Channel()->Lh()->getChannelScaleParam();
  fixedParams.push_back(fixedScaleParam);
  Info << "Fixing scaling parameter " << fixedScaleParam << endmsg;

  const std::vector<std::string> parNames=upar->ParamNames();

  std::vector<std::string>::const_iterator itFix;
  for (itFix=fixedParams.begin(); itFix!=fixedParams.end(); ++itFix){
    //check if name exisists
    if(std::find(parNames.begin(), parNames.end(), (*itFix)) != parNames.end()) upar->Fix( (*itFix) );
    else{
      if(mustMatch){
	Alert << "parameter with name\t" << (*itFix) <<"\tdoes not exist!!!" << endmsg;
	exit(0);
      }
      else Warning << "parameter with name\t" << (*itFix) <<"\tdoes not exist!!!" << endmsg; 
    }
  }
}

void AppBase::fixAllReleaseScaleParams(std::shared_ptr<AbsPawianParameters> upar){
  //  std::string scaleParam = GlobalEnv::instance()->Channel()->Lh()->getChannelScaleParam() + "Other";
  std::string scaleParam = GlobalEnv::instance()->Channel()->Lh()->getChannelScaleParam(); 
  const std::vector<std::string> parNames=upar->ParamNames();

  std::vector<std::string>::const_iterator itFix;
  for (itFix=parNames.begin(); itFix!=parNames.end(); ++itFix){
    if (scaleParam != (*itFix)) upar->Fix( (*itFix) );
    else upar->Release((*itFix)  );
  }

}

bool AppBase::calcAndSendClientLh(NetworkClient& theClient, std::shared_ptr<AbsPawianParameters> startParams, ChannelID channelID){

  std::shared_ptr<AbsPawianParameters> currentFitParams= std::shared_ptr<AbsPawianParameters>(startParams->Clone());
  
  while(true){
    if(!theClient.WaitForParams()) return false;

    const std::vector<double> currentParamVec=theClient.GetParams();
    currentFitParams->SetAllValues(currentParamVec);

    LHData theLHData;
    GlobalEnv::instance()->Channel(channelID)->Lh()->calcLogLhDataClient(currentFitParams, theLHData);

    if(!theClient.SendLH(theLHData.logLH_data, theLHData.LH_mc)) return false;
  }
  return true;
}

void AppBase::fitServerMode(std::shared_ptr<AbsPawianParameters> upar){
  double evtWeightSumData=0;
  ChannelEnvList channelEnvs=GlobalEnv::instance()->ChannelEnvs();
  std::map<short, std::tuple<long, double, long> > numEventMap;

  for(auto it=channelEnvs.begin();it!=channelEnvs.end();++it){
    const std::string datFile=(*it).first->parser()->dataFile();
    const std::string mcFile=(*it).first->parser()->mcFile();
    Info << "data file: " << datFile ;  // << endmsg;
    Info << "mc file: " << mcFile ;  // << endmsg;
    
    int noOfDataEvents =(*it).first->parser()->noOfDataEvts();
    int ratioMcToData=(*it).first->parser()->ratioMcToData();
    
    std::vector<std::string> dataFileNames;
    dataFileNames.push_back(datFile);
    
    std::vector<std::string> mcFileNames;
    mcFileNames.push_back(mcFile);
    
    EventList eventsData;
    readEvents(eventsData, dataFileNames, (*it).first->channelID(), (*it).first->useDataEvtWeight(), 0, noOfDataEvents);
    
    EventList mcData;
    int maxMcEvts=eventsData.size()*ratioMcToData;
    readEvents(mcData, mcFileNames, (*it).first->channelID(), (*it).first->useMCEvtWeight(), 0, maxMcEvts-1);
    
    std::shared_ptr<EvtWeightList> evtWeightListPtr(new EvtWeightList((*it).first->channelID()));
    evtWeightListPtr->read(eventsData, mcData);
    evtWeightSumData+=evtWeightListPtr->NoOfWeightedDataEvts();
    
    numEventMap[(*it).first->channelID()] = std::tuple<long, double,long>(eventsData.size(), evtWeightListPtr->NoOfWeightedDataEvts(), mcData.size());
    }

  std::shared_ptr<AbsFcn> absFcn;
  std::shared_ptr<NetworkServer> theServer(new NetworkServer(GlobalEnv::instance()->parser()->serverPort(), GlobalEnv::instance()->parser()->noOfClients(), numEventMap, GlobalEnv::instance()->parser()->clientNumberWeights()));
  //  theServer->WaitForFirstClientLogin();

  absFcn=std::shared_ptr<AbsFcn>(new PwaFcnServer(theServer));
  theServer->WaitForFirstClientLogin();

  std::shared_ptr<AbsPawianMinimizer> absMinimizerPtr;
  if(GlobalEnv::instance()->parser()->mode()=="server") absMinimizerPtr=std::shared_ptr<AbsPawianMinimizer>(new MinuitMinimizer(absFcn, upar));
  else if (GlobalEnv::instance()->parser()->mode()=="evoserver") absMinimizerPtr=std::shared_ptr<AbsPawianMinimizer>(new EvoMinimizer(absFcn, upar, GlobalEnv::instance()->parser()->evoPopulation(), GlobalEnv::instance()->parser()->evoIterations()));
  else{
    Alert << "only the options server or evoserver are supported for the fitServerMode" << endmsg;
    Alert << "thus " << GlobalEnv::instance()->parser()->mode() << " is not supported" << endmsg;
    exit(1);  
  }
  
  absMinimizerPtr->minimize();
  
  absMinimizerPtr->printFitResult(evtWeightSumData);
  absMinimizerPtr->dumpFitResult();

  theServer->BroadcastClosingMessage();
  Info << "Closing server." << endmsg;
}


void AppBase::fitNonServerMode(std::shared_ptr<AbsPawianParameters> upar, double evtWeightSumData){
  std::shared_ptr<AbsFcn> absFcn(new PwaFcnBase());
  std::shared_ptr<AbsPawianMinimizer> absMinimizerPtr;
  if(GlobalEnv::instance()->parser()->mode()=="pwa") absMinimizerPtr=std::shared_ptr<AbsPawianMinimizer>(new MinuitMinimizer(absFcn, upar));
  else if (GlobalEnv::instance()->parser()->mode()=="evo") absMinimizerPtr=std::shared_ptr<AbsPawianMinimizer>(new EvoMinimizer(absFcn, upar, GlobalEnv::instance()->parser()->evoPopulation(), GlobalEnv::instance()->parser()->evoIterations()));
  else{
    Alert << "fitNonServerMode only the options pwa or evo are supported for the fitNonServerMode" << endmsg;
    Alert << "thus " << GlobalEnv::instance()->parser()->mode() << " is not supported" << endmsg;
    exit(1);  
  }
  
  absMinimizerPtr->minimize();
  absMinimizerPtr->printFitResult(evtWeightSumData);
  absMinimizerPtr->dumpFitResult();
}

void AppBase::fitClientMode(std::shared_ptr<AbsPawianParameters> theStartparams){
  //  ChannelEnvList channelEnvs=GlobalEnv::instance()->ChannelEnvs();
  std::ostringstream portStringStream;
  portStringStream << GlobalEnv::instance()->parser()->serverPort();
    NetworkClient theClient(GlobalEnv::instance()->parser()->serverAddress(), portStringStream.str());
  if(!theClient.Login()){
    Alert << "login of the client failed!!!" << endmsg;
    exit(1);
  }

  ChannelID channelID = theClient.channelID();
  bool cacheAmps = GlobalEnv::instance()->Channel(channelID)->parser()->cacheAmps();
  Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
  if (cacheAmps) GlobalEnv::instance()->Channel(channelID)->Lh()->cacheAmplitudes();

  const std::string datFile=GlobalEnv::instance()->Channel(channelID)->parser()->dataFile();
  const std::string mcFile=GlobalEnv::instance()->Channel(channelID)->parser()->mcFile();
  Info << "data file: " << datFile ;  // << endmsg;
  Info << "mc file: " << mcFile ;  // << endmsg;
  
  std::vector<std::string> dataFileNames;
  dataFileNames.push_back(datFile);

  std::vector<std::string> mcFileNames;
  mcFileNames.push_back(mcFile);

  EventList eventsDataClient;
  readEvents(eventsDataClient, dataFileNames, channelID, GlobalEnv::instance()->Channel(channelID)->useDataEvtWeight(), theClient.GetEventLimits()[0], theClient.GetEventLimits()[1]);

  EventList mcDataClient;
  readEvents(mcDataClient, mcFileNames, channelID, GlobalEnv::instance()->Channel(channelID)->useMCEvtWeight(), theClient.GetEventLimits()[2], theClient.GetEventLimits()[3]);

  std::shared_ptr<EvtDataBaseList> eventListPtr(new EvtDataBaseList(channelID));
  eventListPtr->read(eventsDataClient, mcDataClient);

  GlobalEnv::instance()->Channel(channelID)->Lh()->setDataVec(eventListPtr->getDataVecs());
  GlobalEnv::instance()->Channel(channelID)->Lh()->setMcVec(eventListPtr->getMcVecs());

  calcAndSendClientLh(theClient, theStartparams, channelID);
}

void AppBase::addChannelEnvs(int argcWoCfgFile, char** argvWoCfgFile){
  int argcWCfgFile=argcWoCfgFile+2;
  char* argvWCfgFile[argcWCfgFile];
  for (int i=0; i<argcWoCfgFile ; ++i){
    if (i==0) argvWCfgFile[i]=argvWoCfgFile[i];
    else argvWCfgFile[i+2]=argvWoCfgFile[i];
  }

  std::vector<std::string> pbarpCfgs = GlobalEnv::instance()->parser()->pbarpCfgs();
  loopChannelEnvFactory(argcWCfgFile, argvWCfgFile, pbarpCfgs, AbsChannelEnv::CHANNEL_PBARP);
  std::vector<std::string> epemCfgs = GlobalEnv::instance()->parser()->epemCfgs();
  loopChannelEnvFactory(argcWCfgFile, argvWCfgFile, epemCfgs, AbsChannelEnv::CHANNEL_EPEM);
}

void AppBase::loopChannelEnvFactory(int argcWCfgFile, char** argvWCfgFile, std::vector<std::string>& reactionCfgs, short channelType){
  for(auto it=reactionCfgs.begin(); it!=reactionCfgs.end();++it){
     argvWCfgFile[1]=(char*)"-c";
     argvWCfgFile[2]=(char*)(*it).c_str();
     for (int i=0; i<argcWCfgFile ; ++i){
       Info << "argvWCfgFile[" << i << "]= " << argvWCfgFile[i] << endmsg;
     }
     std::shared_ptr<AbsChannelEnv> channelEnv;
     if(channelType==AbsChannelEnv::CHANNEL_PBARP){
       pbarpParser* currentParser = new pbarpParser(argcWCfgFile, argvWCfgFile);
       channelEnv = std::shared_ptr<AbsChannelEnv>(new PbarpChannelEnv(currentParser));
     }
     else if(channelType==AbsChannelEnv::CHANNEL_EPEM){
        epemParser* currentParser = new epemParser(argcWCfgFile, argvWCfgFile);
	channelEnv = std::shared_ptr<AbsChannelEnv>(new EpemChannelEnv(currentParser));
   }
     GlobalEnv::instance()->AddEnv(channelEnv, channelType);
  }
}
