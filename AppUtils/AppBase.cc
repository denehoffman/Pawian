//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universit??t Bochum 				  //
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
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/PwaGen.hh"
#include "PwaUtils/AbsHist.hh"
#include "qaErrorExtract/WaveContribution.hh"
#include "qaErrorExtract/ScatteringContribution.hh"
#include "PwaUtils/NetworkServer.hh"
#include "PwaUtils/NetworkClient.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/EvtDataListFactory.hh"

#include "FitParams/AbsPawianParamStreamer.hh"
#include "FitParams/PwaCovMatrix.hh"

#include "ConfigParser/ParserBase.hh"
#include "ConfigParser/pbarpParser.hh"
#include "ConfigParser/epemParser.hh"
#include "ConfigParser/resParser.hh"
#include "ConfigParser/ggParser.hh"
#include "ConfigParser/pipiScatteringParser.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Event/Event.hh"
#include "Event/EventReaderDefault.hh"
#include "Event/EventReaderScattering.hh"

#include "MinFunctions/PwaFcnBase.hh"
#include "MinFunctions/PwaFcnServer.hh"
#include "MinFunctions/AbsPawianMinimizer.hh"
#include "MinFunctions/EvoMinimizer.hh"
#include "MinFunctions/MinuitMinimizer.hh"

#include "FitParams/AbsPawianParameters.hh"
#include "FitParams/ParamDepHandler.hh"

#include "pbarpUtils/PbarpChannelEnv.hh"
#include "epemUtils/EpemChannelEnv.hh"
#include "resUtils/ResChannelEnv.hh"
#include "ggUtils/GGChannelEnv.hh"
#include "pipiScatteringUtils/PiPiScatteringChannelEnv.hh"

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

void AppBase::dumpRandomParams(){
  std::stringstream randomparamsname;
  randomparamsname << "randomParams" << GlobalEnv::instance()->outputFileNameSuffix() << ".dat";
  std::ofstream theStreamRandom ( randomparamsname.str().c_str() );
  
  std::shared_ptr<AbsPawianParameters> randomParams=GlobalEnv::instance()->randomPawianParams();
  randomParams->print(theStreamRandom);
}

void AppBase::generate(std::shared_ptr<AbsPawianParameters> theParams){
    std::shared_ptr<PwaGen> pwaGenPtr(new PwaGen());
    theParams->print(std::cout);
    pwaGenPtr->generate(GlobalEnv::instance()->Channel()->Lh(), theParams);
}

void AppBase::readEvents(EventList& theEventList, std::vector<std::string>& fileNames, 
			 ChannelID channelID, bool withEvtWeight, int evtStart, int evtStop){
  std::vector<std::string> completeFileNames;
  std::vector<std::string>::iterator itStr;
  for(itStr=fileNames.begin(); itStr!=fileNames.end(); ++itStr){
    std::string currentFile=GlobalEnv::instance()->Channel(channelID)->parser()->prePathDataFiles()+(*itStr);
    completeFileNames.push_back(currentFile);
  }
  int noFinalStateParticles=GlobalEnv::instance()->Channel(channelID)->noFinalStateParticles();
  std::vector< std::shared_ptr<MassRangeCut> > massRangeCuts =
    GlobalEnv::instance()->Channel(channelID)->massRangeCuts();
  EventReaderDefault eventReader(completeFileNames, noFinalStateParticles, 0, withEvtWeight);
  eventReader.setUnit(GlobalEnv::instance()->Channel(channelID)->parser()->unitInFile());
  eventReader.setOrder(GlobalEnv::instance()->Channel(channelID)->parser()->orderInFile());

  if(GlobalEnv::instance()->Channel(channelID)->useMassRange()){
   eventReader.setMassRange(massRangeCuts);
  }

  eventReader.fill(theEventList, evtStart, evtStop);

  unsigned int numOfEvtsInFile=theEventList.size();
  InfoMsg << "File has " << numOfEvtsInFile << " events." << endmsg;
  if (numOfEvtsInFile==0){
    Alert << "0 events in the file; abort!!!" << endmsg;
    exit(1);
  }

  Event* currentNextEvent = theEventList.nextEvent();
  if (0!=currentNextEvent) { 
    InfoMsg << " Each event has "
	    <<  theEventList.nextEvent()->size() << " final state particles.\n" << endmsg;
  } else {
    WarningMsg << "\n currentNextEvent does not exist.\n" << endmsg;
  }

  theEventList.rewind();

  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = theEventList.nextEvent()) != 0 && evtCount < 10) {
    InfoMsg << endmsg;
    for(int i=0; i<noFinalStateParticles; ++i){
      InfoMsg << (*anEvent->p4(i)) << "\tm = " << anEvent->p4(i)->Mass() << endmsg;
    }
    ++evtCount;
  }
  InfoMsg << endmsg;
  theEventList.rewind();
}

void AppBase::readScatteringEvents(EventList& theEventList, std::vector<std::string>& fileNames, 
				   ChannelID channelID, int evtStart, int evtStop){
  std::vector<std::string> completeFileNames;
  std::vector<std::string>::iterator itStr;
  for(itStr=fileNames.begin(); itStr!=fileNames.end(); ++itStr){
    std::string currentFile=GlobalEnv::instance()->Channel(channelID)->parser()->prePathDataFiles()+(*itStr);
    completeFileNames.push_back(currentFile);
  }

  std::vector< std::shared_ptr<MassRangeCut> > massRangeCuts =
    GlobalEnv::instance()->Channel(channelID)->massRangeCuts();
  EventReaderScattering evtScatterReader(completeFileNames, 2, 0, false);
  if(GlobalEnv::instance()->Channel(channelID)->useMassRange()){
    evtScatterReader.setMassRange(massRangeCuts);
  }

  evtScatterReader.fill(theEventList, evtStart, evtStop);
}

void AppBase::createLhObjects(){

  ChannelEnvList channelEnvs=GlobalEnv::instance()->ChannelEnvs();
  for(auto it=channelEnvs.begin();it!=channelEnvs.end();++it){
    std::string prodFormalism=(*it).first->parser()->productionFormalism();
    InfoMsg <<"create likelihood object for channelType " << (*it).first->channelType()
	    <<"\t and channelID: " << (*it).first->channelID() << endmsg;
    std::shared_ptr<AbsLh> theLhPtr=LhFactory::instance()->getLh((*it).first->channelType(),
								 (*it).first->channelID(), prodFormalism );
    (*it).first->SetLh(theLhPtr);
  }
  
  // Generate the full parameter set using the likelihood list
  //  GlobalEnv::instance()->CreateDefaultParameterSet();
}

void AppBase::qaMode(std::shared_ptr<AbsPawianParameters> startParams, double evtWeightSumData){

  std::ostringstream qaSummaryFileName;
  std::string outputFileNameSuffix= GlobalEnv::instance()->outputFileNameSuffix();
  qaSummaryFileName << "qaSummary" << outputFileNameSuffix << ".dat";
  std::ofstream theQaStream ( qaSummaryFileName.str().c_str() );

  int noOfFreeFitParams=startParams->VariableParameters();
  if(evtWeightSumData<=(noOfFreeFitParams+1)) {
    WarningMsg << "number of data events less or equal to the number of free parameters!!!"
	       << "\n be careful with the fit result!!!"
	       << endmsg; 
  }  
  
  double theLh=GlobalEnv::instance()->Channel()->Lh()->calcLogLh(startParams);
  double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
  double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
  double AICccDenom = evtWeightSumData-noOfFreeFitParams-1;
  if ( fabs(AICccDenom) < 1.e-10) AICccDenom=1.e-10;
  double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/AICccDenom;

  std::shared_ptr<PwaCovMatrix> thePwaCovMatrix(new PwaCovMatrix);
  if(GlobalEnv::instance()->parser()->calcContributionError()){
    std::string serializationFileName = GlobalEnv::instance()->serializationFileName();
    std::ifstream serializationStream(serializationFileName.c_str());
    
    if(!serializationStream.is_open()){
      Alert << "Could not open serialization file." << endmsg;
      exit(0);
    }
    boost::archive::text_iarchive boostInputArchive(serializationStream);
    boostInputArchive >> *thePwaCovMatrix;
  }

  if(GlobalEnv::instance()->Channel()->channelType() == AbsChannelEnv::CHANNEL_PIPISCATTERING){
    if(GlobalEnv::instance()->parser()->calcContributionError()){
      std::shared_ptr<ScatteringContribution> theScatteringContribution(
	  new ScatteringContribution(GlobalEnv::instance()->Channel()->Lh(), startParams, thePwaCovMatrix));
      theScatteringContribution->CalcError();
    }

    std::shared_ptr<AbsHist> histPtr1 = GlobalEnv::instance()->Channel()->CreateHistInstance();
    histPtr1->fillFromLhData(GlobalEnv::instance()->Channel()->Lh(), startParams);
    InfoMsg << "chi2\t" << theLh << endmsg;
    theQaStream << "chi2\t" << theLh << "\n";
    return;
  }
  

  std::shared_ptr<WaveContribution> theWaveContribution;
  if(GlobalEnv::instance()->parser()->calcContributionError()){
    theWaveContribution = std::shared_ptr<WaveContribution>
      (new WaveContribution(GlobalEnv::instance()->Channel()->Lh(), startParams, thePwaCovMatrix));
  }
  else{
    theWaveContribution = std::shared_ptr<WaveContribution>
      (new WaveContribution(GlobalEnv::instance()->Channel()->Lh(), startParams));
  }
  std::pair<double, double> contValue = theWaveContribution->CalcContribution();
  std::vector<std::pair<std::string,std::pair<double,double>>> singleContValues 
    = theWaveContribution->CalcSingleContributions();

  InfoMsg << "logLh\t" << theLh << endmsg;
  theQaStream << "logLh\t" << theLh << "\n";

  InfoMsg << "noOfFreeFitParams:\t" << noOfFreeFitParams;
  theQaStream << "noOfFreeFitParams:\t" << noOfFreeFitParams << "\n";

  InfoMsg << "BIC:\t" << BICcriterion << endmsg;
  theQaStream << "BIC:\t" << BICcriterion << "\n";

  InfoMsg << "AICa:\t" << AICcriterion << endmsg;
  theQaStream << "AICa:\t" << AICcriterion << "\n";

  InfoMsg << "AICc:\t" << AICccriterion << endmsg;
  theQaStream << "AICc:\t" << AICccriterion << "\n";

  InfoMsg << "No of data events without weight " 
	  << GlobalEnv::instance()->Channel()->Lh()->getDataVec().size() << endmsg;
  theQaStream << "No of data events without weight " 
	      << GlobalEnv::instance()->Channel()->Lh()->getDataVec().size() << "\n";

  InfoMsg << "No of data events with weight " << evtWeightSumData << endmsg;
  theQaStream << "No of data events with weight " << evtWeightSumData << "\n";

  InfoMsg << "No of MC events " << GlobalEnv::instance()->Channel()->Lh()->getMcVec().size() << endmsg;
  theQaStream << "No of MC events " << GlobalEnv::instance()->Channel()->Lh()->getMcVec().size() << "\n";

  double scaleFactor = evtWeightSumData/GlobalEnv::instance()->Channel()->Lh()->getMcVec().size();
  InfoMsg << "scaling factor " << scaleFactor << endmsg;
  theQaStream << "scaling factor " << scaleFactor << "\n";

  InfoMsg << "no of fitted events with scaling factor: " << contValue.first*scaleFactor << endmsg;
  theQaStream << "no of fitted events with scaling factor: " << contValue.first*scaleFactor << "\n";

  InfoMsg << "Selected wave contribution:\t" << contValue.first << " +- " << contValue.second;
  theQaStream << "Selected wave contribution:\t" << contValue.first << " +- " << contValue.second <<  "\n";

  std::vector<std::pair<std::string,std::pair<double,double>>>::iterator it;
  for(it=singleContValues.begin(); it!=singleContValues.end(); ++it) {
    InfoMsg << "Single wave contribution " << (*it).first << "\t" << (*it).second.first 
	    << " +- " << (*it).second.second;
    theQaStream << "Single wave contribution " << (*it).first << "\t" << (*it).second.first 
		<< " +- " << (*it).second.second <<  "\n";
  }

  InfoMsg << endmsg;
  theQaStream.close();

  std::shared_ptr<AbsHist> histPtr = GlobalEnv::instance()->Channel()->CreateHistInstance();
  histPtr->fillFromLhData(GlobalEnv::instance()->Channel()->Lh(), startParams);
}

void AppBase::qaModeSimple(EventList& dataEventList, EventList& mcEventList, 
			   std::shared_ptr<AbsPawianParameters> startParams){
  std::shared_ptr<EvtDataBaseList> evtDataBaseList =
    EvtDataListFactory::instance()->evtDataListPtr(GlobalEnv::instance()->Channel());

  int noOfFreeFitParams=startParams->VariableParameters();
  std::shared_ptr<AbsLh> absLh=GlobalEnv::instance()->Channel()->Lh();
  LHData theLHData;
  std::shared_ptr<WaveContribution> theWaveContribution(
      new WaveContribution(GlobalEnv::instance()->Channel()->Lh(), startParams));

  std::shared_ptr<AbsPawianParameters> currentParams = 
    std::shared_ptr<AbsPawianParameters>(startParams->Clone());

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
    int dataPoint=1;
    while ((anEvent = dataEventList.nextEvent())){
      EvtData* currentDataEvt=evtDataBaseList->convertEvent(anEvent, evtCount);
      absLh->addDataToLogLh(currentDataEvt, currentParams, theLHData);
      histPtr->fillEvt(currentDataEvt, currentDataEvt->evtWeight, "data", dataPoint);
      evtWeightSumData += currentDataEvt->evtWeight;
      delete currentDataEvt;
      evtCount++;
      dataPoint++;
      if (evtCount%1000 == 0) {
	InfoMsg << evtCount << " data events calculated" << endmsg;
      }
    }

    //loop over mc events
    int evtCountMc = 0;

    double integralFitWeight=0.;

    mcEventList.rewind();
    dataPoint=1;
    while ((anEvent = mcEventList.nextEvent())){
      EvtData* currentMcEvt=evtDataBaseList->convertEvent(anEvent, evtCount);
      double currentIntensity=absLh->addMcToLogLh(currentMcEvt,currentParams, theLHData);
      histPtr->fillEvt(currentMcEvt, 1., "mc", dataPoint);
      histPtr->fillEvt(currentMcEvt, currentIntensity, "fit", dataPoint);

      integralFitWeight+=currentIntensity;

      delete currentMcEvt;
      evtCount++;
      evtCountMc++;
      dataPoint++;
      if (evtCountMc%1000 == 0) {
	InfoMsg << evtCountMc << " MC events calculated" << endmsg;
      }
    }

    double scaleFactor=theLHData.weightSum/theLHData.num_mc;
    histPtr->scaleFitHists(scaleFactor);

    if(i==-1){
      double theLh=absLh->mergeLogLhData(theLHData);
      double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
      double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
      double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1) /
	(evtWeightSumData-noOfFreeFitParams-1);
      
      InfoMsg << "logLh\t" << theLh << endmsg;
      theQaStream << "logLh\t" << theLh << "\n";
      
      InfoMsg << "noOfFreeFitParams:\t" << noOfFreeFitParams;
      theQaStream << "noOfFreeFitParams\t" << noOfFreeFitParams << "\n";
      
      InfoMsg << "BIC:\t" << BICcriterion << endmsg;
      theQaStream << "BIC:\t" << BICcriterion << "\n";
      
      InfoMsg << "AICa:\t" << AICcriterion << endmsg;
      theQaStream << "AICa:\t" << AICcriterion << "\n";
      
      InfoMsg << "AICc:\t" << AICccriterion << endmsg;
      theQaStream << "AICc:\t" << AICccriterion << "\n";
      
      InfoMsg << "No of data events without weight " << integralDataWoWeight << endmsg;
      theQaStream << "No of data events without weight " << integralDataWoWeight << "\n";
      
      InfoMsg << "No of data events with weight " << evtWeightSumData << endmsg;
      theQaStream << "No of data events with weight " << evtWeightSumData << "\n";
      
      InfoMsg << "No of MC events " << theLHData.num_mc << endmsg;
      theQaStream << "No of MC events " << theLHData.num_mc << "\n";
      
      InfoMsg << "scaling factor " << scaleFactor << endmsg;
      theQaStream << "scaling factor " << scaleFactor << "\n";
      
      InfoMsg << "no of fitted events with scaling factor: " << integralFitWeight*scaleFactor;
      theQaStream << "no of fitted events with scaling factor: " << integralFitWeight*scaleFactor << "\n";
    }
    else{ //i>-1
      InfoMsg << "contribution no " << i << endmsg;
      InfoMsg << "contribution name " << contributionName << endmsg;
      InfoMsg << "No of data events without weight " << integralDataWoWeight << endmsg;
      InfoMsg << "No of data events with weight " << evtWeightSumData << endmsg;
      InfoMsg << "No of MC events " << theLHData.num_mc << endmsg;
      InfoMsg << "no of fitted events with scaling factor: " << integralFitWeight*scaleFactor << endmsg;
      theQaStream << "contribution no " << i << "\n";
      theQaStream << "contribution name " << contributionName << "\n";
      theQaStream << "no of fitted events with scaling factor: " << integralFitWeight*scaleFactor << "\n";
    }
  } // loop over contributions
  InfoMsg << endmsg;
  theQaStream.close();
}

void AppBase::qaModeEffCorrection(EventList& dataEventList, EventList& mcEventList,
				  EventList& truthEventList, 
				  std::shared_ptr<AbsPawianParameters> startParams) {
  std::shared_ptr<EvtDataBaseList> evtDataBaseList(new EvtDataBaseList(0));
  std::shared_ptr<AbsLh> absLh=GlobalEnv::instance()->Channel()->Lh();
  LHData theLHData;
  std::shared_ptr<AbsPawianParameters> currentParams = 
    std::shared_ptr<AbsPawianParameters>(startParams->Clone());

  std::shared_ptr<AbsHist> histPtr= GlobalEnv::instance()->Channel()->CreateHistInstance("", true);

  absLh->updateFitParams(currentParams);
  
  //loop over data events
  Event* anEvent;
  int evtCount = 0;
  
  double evtWeightSumData=0.;
  
  dataEventList.rewind();
  int dataPoint=1;
  while ((anEvent = dataEventList.nextEvent())){
    EvtData* currentDataEvt=evtDataBaseList->convertEvent(anEvent, evtCount);
    absLh->addDataToLogLh(currentDataEvt, currentParams, theLHData);
    histPtr->fillEvt(currentDataEvt, currentDataEvt->evtWeight, "data", dataPoint);
    evtWeightSumData += currentDataEvt->evtWeight;
    delete currentDataEvt;
    evtCount++;
    dataPoint++;
    if (evtCount%1000 == 0) {
      InfoMsg << evtCount << " data events calculated" << endmsg;
    }
  }
  
  //loop over mc events
  int evtCountMc = 0;
  
  double integralFitWeight=0.;
  
  mcEventList.rewind();
  dataPoint=1;
  while ((anEvent = mcEventList.nextEvent())) {
    EvtData* currentMcEvt=evtDataBaseList->convertEvent(anEvent, evtCount);
    double currentIntensity=absLh->addMcToLogLh(currentMcEvt,currentParams, theLHData);
    histPtr->fillEvt(currentMcEvt, 1., "mc", dataPoint);
    histPtr->fillEvt(currentMcEvt, currentIntensity, "fit", dataPoint);
    
    integralFitWeight+=currentIntensity;
    
    delete currentMcEvt;
    evtCount++;
    evtCountMc++;
    dataPoint++;
    if (evtCountMc%1000 == 0) {
      InfoMsg << evtCountMc << " MC events calculated" << endmsg;
    }
  }
  
  //loop over truth events
  int evtCountTruth = 0;
  
  double integralTruthFitWeight=0.;
  
  truthEventList.rewind();
  dataPoint=1;
  while ((anEvent = truthEventList.nextEvent())) {
    EvtData* currentTruthEvt=evtDataBaseList->convertEvent(anEvent, evtCount);
    double currentIntensity=absLh->calcEvtIntensity(currentTruthEvt, currentParams);
    histPtr->fillEvt(currentTruthEvt, 1., "truthWoWeight", dataPoint);
    histPtr->fillEvt(currentTruthEvt, currentIntensity, "truthWWeight", dataPoint);
    
    integralTruthFitWeight+=currentIntensity;
    
    delete currentTruthEvt;
    evtCount++;
    evtCountTruth++;
    dataPoint++;
    if (evtCountTruth%1000 == 0) {
      InfoMsg << evtCountTruth << " Truth events calculated" << endmsg ;
    }
  }
  
  double scaleFactor=theLHData.weightSum/theLHData.num_mc;
  histPtr->scaleFitHists(scaleFactor);
}


void AppBase::plotMode(EventList& dataEventList, EventList& mcEventList, 
		       std::shared_ptr<EvtDataBaseList> evtDataBaseList){

  std::shared_ptr<AbsHist> histPtr = GlobalEnv::instance()->Channel()->CreateHistInstance();

  //loop over data events
  Event* anEvent;
  int evtCount = 0;
  double evtWeightSumData=0.;
  dataEventList.rewind();
  int dataPoint=1;
  while ((anEvent = dataEventList.nextEvent())){
    EvtData* currentDataEvt=evtDataBaseList->convertEvent(anEvent, anEvent->eventNo());
    histPtr->fillEvt(currentDataEvt, currentDataEvt->evtWeight, "data", dataPoint);
    if(GlobalEnv::instance()->Channel()->channelType() == AbsChannelEnv::CHANNEL_PIPISCATTERING){
      histPtr->fillEvt(currentDataEvt, currentDataEvt->evtWeight, "fit", dataPoint);
    }
    evtWeightSumData+=currentDataEvt->evtWeight;
    delete currentDataEvt;
    evtCount++;
    dataPoint++;
    if (evtCount%1000 == 0) {
      InfoMsg << evtCount << " data events calculated" << endmsg;
    }
  }

  //loop over mc events
  int evtCountMc = 0;
  if(GlobalEnv::instance()->Channel()->channelType() != AbsChannelEnv::CHANNEL_PIPISCATTERING){
    mcEventList.rewind();
    dataPoint=1;
    while ((anEvent = mcEventList.nextEvent())){
      // EvtData* currentMcEvt=evtDataBaseList->convertEvent(anEvent, evtCount);
      EvtData* currentMcEvt=evtDataBaseList->convertEvent(anEvent);
      histPtr->fillEvt(currentMcEvt, 1., "fit", dataPoint);
      histPtr->fillEvt(currentMcEvt, 1., "mc", dataPoint);
      delete currentMcEvt;
      evtCount++;
      evtCountMc++;
      dataPoint++;
      if (evtCountMc%1000 == 0) {
	InfoMsg << evtCountMc << " MC events calculated" << endmsg;
      }
    }
  }

  double histScaleFactor=evtWeightSumData/evtCountMc;
  histPtr->scaleFitHists(histScaleFactor);
}

std::shared_ptr<AbsPawianParameters> AppBase::streamPawianParams(){
  std::string paramStreamerPath=GlobalEnv::instance()->parser()->fitParamFile();
  AbsPawianParamStreamer thePawianStreamer(paramStreamerPath);
  return thePawianStreamer.paramList();
}

void AppBase::fixParams(std::shared_ptr<AbsPawianParameters> upar, 
			std::vector<std::string> fixedParams, bool mustMatch){

  // Evaluate parameter dependencies and add fixes
  std::vector<std::string> paramDependencies = GlobalEnv::instance()->paramDependencies();

  ParamDepHandler::instance()->Fill(paramDependencies, upar);
  std::vector<std::string> dependentParameters = ParamDepHandler::instance()->DependentParameterNames();

  std::vector<std::string>::const_iterator strIt;
  for(strIt=dependentParameters.begin(); strIt!=dependentParameters.end(); ++strIt){
    std::vector<unsigned int> depIdVec=ParamDepHandler::instance()->RefIds(*strIt);
    upar->SetDependencies( *strIt, depIdVec);
  }

  // Always fix the primary channel's scaling parameters
  std::string fixedScaleParam = GlobalEnv::instance()->Channel()->Lh()->getChannelScaleParam();
  fixedParams.push_back(fixedScaleParam);
  InfoMsg << "Fixing scaling parameter " << fixedScaleParam << endmsg;

  const std::vector<std::string> parNames=upar->ParamNames();

  std::vector<std::string>::const_iterator itFix;
  for (itFix=fixedParams.begin(); itFix!=fixedParams.end(); ++itFix){
    //check if name exisists
    if(std::find(parNames.begin(), parNames.end(), (*itFix)) != parNames.end()) upar->Fix( (*itFix) );
    else{
      if(mustMatch) {
	Alert << "parameter with name\t" << (*itFix) <<"\tdoes not exist!!!" << endmsg;
	exit(0);
      } else {
	WarningMsg << "parameter with name\t" << (*itFix) <<"\tdoes not exist!!!" << endmsg; 
      }
    }
  }

  if(GlobalEnv::instance()->parser()->fixAllPhases()){
    for (itFix=parNames.begin(); itFix!=parNames.end(); ++itFix){
      if( itFix->substr( itFix->length() - 3 ) == "Phi"){
	upar->Fix( (*itFix) );
	InfoMsg << "Fixing phase " << (*itFix) << endmsg;
      }
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

bool AppBase::calcAndSendClientLh(NetworkClient& theClient, 
				  std::shared_ptr<AbsPawianParameters> startParams, 
				  ChannelID channelID){

  std::shared_ptr<AbsPawianParameters> currentFitParams = 
    std::shared_ptr<AbsPawianParameters>(startParams->Clone());
  
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
    InfoMsg << "data file: " << datFile << endmsg;
    InfoMsg << "mc file: " << mcFile << endmsg;
    
    int noOfDataEvents =(*it).first->parser()->noOfDataEvts();
    int ratioMcToData=(*it).first->parser()->ratioMcToData();
    
    std::vector<std::string> dataFileNames;
    dataFileNames.push_back(datFile);
    
    std::vector<std::string> mcFileNames;
    mcFileNames.push_back(mcFile);

    double noOfMcEvts=0.;

    EventList eventsData;
    EventList mcData;
    if(GlobalEnv::instance()->Channel((*it).first->channelID())->channelType() == 
       AbsChannelEnv::CHANNEL_PIPISCATTERING){
      readScatteringEvents(eventsData, dataFileNames, (*it).first->channelID(), 0, noOfDataEvents);
    }
    else{
      readEvents(eventsData, dataFileNames, (*it).first->channelID(), 
		 (*it).first->useDataEvtWeight(), 0, noOfDataEvents);
    }

    int noDataEvts=eventsData.size();
    double noOfWeightedDataEvts=EvtDataBaseList::noOfWeightedEvts(eventsData, 
								  (*it).first->channelID(), 
								  noOfDataEvents, 0);
    evtWeightSumData+=noOfWeightedDataEvts; 
    eventsData.removeAndDeleteEvents(0, eventsData.size()-1);        

    if(GlobalEnv::instance()->Channel((*it).first->channelID())->channelType() 
       != AbsChannelEnv::CHANNEL_PIPISCATTERING){
      int maxMcEvts=noDataEvts*ratioMcToData;
      
      readEvents(mcData, mcFileNames, (*it).first->channelID(), (*it).first->useMCEvtWeight(), 0, maxMcEvts-1);
      noOfMcEvts=EvtDataBaseList::noOfWeightedEvts(mcData, (*it).first->channelID(), maxMcEvts-1, 0);    
    }


    numEventMap[(*it).first->channelID()] = 
      std::tuple<long, double,long>(noDataEvts, noOfWeightedDataEvts, mcData.size());
    mcData.removeAndDeleteEvents(0, mcData.size()-1);

    if(noOfWeightedDataEvts<5.){
      Alert << "number of weighted data events too small: " << noOfWeightedDataEvts << endmsg;
      exit(1);
    }
    if(GlobalEnv::instance()->Channel((*it).first->channelID())->channelType() 
       != AbsChannelEnv::CHANNEL_PIPISCATTERING && noOfMcEvts<10.){
      if(noOfMcEvts<10.){
	Alert << "number of weighted Monte Carlo events too small: " << noOfMcEvts << endmsg;
	exit(1);
      }
    }
  }
 
  std::shared_ptr<AbsFcn> absFcn;
  std::shared_ptr<NetworkServer> theServer(new NetworkServer(GlobalEnv::instance()->parser()->serverPort(), 
							     GlobalEnv::instance()->parser()->noOfClients(), 
							     numEventMap, 
							     GlobalEnv::instance()->parser()->
							       clientNumberWeights()));
  absFcn=std::shared_ptr<AbsFcn>(new PwaFcnServer(theServer));
  theServer->WaitForFirstClientLogin();

  std::shared_ptr<AbsPawianMinimizer> absMinimizerPtr;
  if(GlobalEnv::instance()->parser()->mode()=="server") 
    absMinimizerPtr = std::shared_ptr<AbsPawianMinimizer>(new MinuitMinimizer(absFcn, upar));
  else if (GlobalEnv::instance()->parser()->mode()=="evoserver") 
    absMinimizerPtr =
      std::shared_ptr<AbsPawianMinimizer>(new EvoMinimizer(absFcn, upar, 
							   GlobalEnv::instance()->parser()->evoPopulation(), 
							   GlobalEnv::instance()->parser()->evoIterations()));
  else{
    Alert << "only the options server or evoserver are supported for the fitServerMode" << endmsg;
    Alert << "thus " << GlobalEnv::instance()->parser()->mode() << " is not supported" << endmsg;
    exit(1);  
  }
  
  absMinimizerPtr->minimize();
  
  absMinimizerPtr->printFitResult(evtWeightSumData);
  absMinimizerPtr->dumpFitResult();

  theServer->BroadcastClosingMessage();
  InfoMsg << "Closing server." << endmsg;
}


void AppBase::fitNonServerMode(std::shared_ptr<AbsPawianParameters> upar, 
			       double evtWeightSumData, double evtWeightSumMc){
    if(evtWeightSumData<10.){
      Alert << "number of wieghted data events too small: " << evtWeightSumData << endmsg;
      exit(1);
    }
    if( (GlobalEnv::instance()->Channel()->channelType() != AbsChannelEnv::CHANNEL_PIPISCATTERING) 
	&& evtWeightSumMc<10.){
      Alert << "number of weighted Monte Carlo events too small: " << evtWeightSumMc << endmsg;
      exit(1);
    }
  std::shared_ptr<AbsFcn> absFcn(new PwaFcnBase());
  std::shared_ptr<AbsPawianMinimizer> absMinimizerPtr;
  if(GlobalEnv::instance()->parser()->mode()=="pwa") 
    absMinimizerPtr=std::shared_ptr<AbsPawianMinimizer>(new MinuitMinimizer(absFcn, upar));
  else if (GlobalEnv::instance()->parser()->mode()=="evo") 
    absMinimizerPtr=std::shared_ptr<AbsPawianMinimizer>(new EvoMinimizer(absFcn, upar, 
	  GlobalEnv::instance()->parser()->evoPopulation(), GlobalEnv::instance()->parser()->evoIterations()));
  else {
    Alert << "fitNonServerMode only the options pwa or evo are supported for the fitNonServerMode" << endmsg;
    Alert << "thus " << GlobalEnv::instance()->parser()->mode() << " is not supported" << endmsg;
    exit(1);  
  }
  
  absMinimizerPtr->minimize();
  absMinimizerPtr->printFitResult(evtWeightSumData);
  absMinimizerPtr->dumpFitResult();
}

void AppBase::fitClientMode(std::shared_ptr<AbsPawianParameters> theStartparams){
  std::ostringstream portStringStream;
  portStringStream << GlobalEnv::instance()->parser()->serverPort();
    NetworkClient theClient(GlobalEnv::instance()->parser()->serverAddress(), portStringStream.str());
  if(!theClient.Login()){
    Alert << "login of the client failed!!!" << endmsg;
    exit(1);
  }

  ChannelID channelID = theClient.channelID();
  bool cacheAmps = GlobalEnv::instance()->Channel(channelID)->parser()->cacheAmps();
  InfoMsg << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
  if (cacheAmps) GlobalEnv::instance()->Channel(channelID)->Lh()->cacheAmplitudes();

  const std::string datFile=GlobalEnv::instance()->Channel(channelID)->parser()->dataFile();
  const std::string mcFile=GlobalEnv::instance()->Channel(channelID)->parser()->mcFile();
  InfoMsg << "data file: " << datFile << endmsg;
  InfoMsg << "mc file: " << mcFile << endmsg;
  
  std::vector<std::string> dataFileNames;
  dataFileNames.push_back(datFile);

  std::vector<std::string> mcFileNames;
  mcFileNames.push_back(mcFile);

  EventList eventsDataClient;
  EventList mcDataClient;
  if(GlobalEnv::instance()->Channel(channelID)->channelType() == AbsChannelEnv::CHANNEL_PIPISCATTERING){
     readScatteringEvents(eventsDataClient, dataFileNames, channelID, theClient.GetEventLimits()[0], 
			  theClient.GetEventLimits()[1]);
   }
  else {
    readEvents(eventsDataClient, dataFileNames, channelID, 
	       GlobalEnv::instance()->Channel(channelID)->useDataEvtWeight(), 
	       theClient.GetEventLimits()[0], theClient.GetEventLimits()[1]);
    
    readEvents(mcDataClient, mcFileNames, channelID, 
	       GlobalEnv::instance()->Channel(channelID)->useMCEvtWeight(), 
	       theClient.GetEventLimits()[2], theClient.GetEventLimits()[3]);
      }

  std::shared_ptr<EvtDataBaseList> eventListPtr
    = EvtDataListFactory::instance()->evtDataListPtr(GlobalEnv::instance()->Channel(channelID));
  eventListPtr->read(eventsDataClient, mcDataClient);

  eventsDataClient.removeAndDeleteEvents(0, eventsDataClient.size()-1);
  mcDataClient.removeAndDeleteEvents(0, mcDataClient.size()-1);

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
  std::vector<std::string> ggCfgs = GlobalEnv::instance()->parser()->ggCfgs();
  loopChannelEnvFactory(argcWCfgFile, argvWCfgFile, ggCfgs, AbsChannelEnv::CHANNEL_GG);
  std::vector<std::string> resCfgs = GlobalEnv::instance()->parser()->resCfgs();
  loopChannelEnvFactory(argcWCfgFile, argvWCfgFile, resCfgs, AbsChannelEnv::CHANNEL_RES);
  std::vector<std::string> pipiScatteringCfgs = GlobalEnv::instance()->parser()->pipiScatteringCfgs();
  loopChannelEnvFactory(argcWCfgFile, argvWCfgFile, pipiScatteringCfgs, AbsChannelEnv::CHANNEL_PIPISCATTERING);
}

void AppBase::loopChannelEnvFactory(int argcWCfgFile, char** argvWCfgFile, 
				    std::vector<std::string>& reactionCfgs, short channelType){
  for(auto it=reactionCfgs.begin(); it!=reactionCfgs.end();++it){
     argvWCfgFile[1]=(char*)"-c";
     argvWCfgFile[2]=(char*)(*it).c_str();
     for (int i=0; i<argcWCfgFile ; ++i){
       InfoMsg << "argvWCfgFile[" << i << "]= " << argvWCfgFile[i] << endmsg;
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
     else if(channelType==AbsChannelEnv::CHANNEL_RES){
        resParser* currentParser = new resParser(argcWCfgFile, argvWCfgFile);
	channelEnv = std::shared_ptr<AbsChannelEnv>(new ResChannelEnv(currentParser));
   }
     else if(channelType==AbsChannelEnv::CHANNEL_GG){
        ggParser* currentParser = new ggParser(argcWCfgFile, argvWCfgFile);
	channelEnv = std::shared_ptr<AbsChannelEnv>(new GGChannelEnv(currentParser));
   }
     else if(channelType==AbsChannelEnv::CHANNEL_PIPISCATTERING){
        pipiScatteringParser* currentParser = new pipiScatteringParser(argcWCfgFile, argvWCfgFile);
	channelEnv = std::shared_ptr<AbsChannelEnv>(new PiPiScatteringChannelEnv(currentParser));
   }
     GlobalEnv::instance()->AddEnv(channelEnv, channelType);
  }
}
