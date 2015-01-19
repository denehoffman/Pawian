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

#include "PwaUtils/AppBase.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/PwaGen.hh"
#include "ConfigParser/ParserBase.hh"
#include "PwaUtils/AbsHist.hh"
#include "PwaUtils/WaveContribution.hh"
#include "PwaUtils/PwaCovMatrix.hh"
#include "PwaUtils/NetworkClient.hh"
#include "PwaUtils/EvtDataBaseList.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Event/Event.hh"
#include "Event/EventReaderDefault.hh"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserCovariance.h"

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

    fitParams paramVal = GlobalEnv::instance()->DefaultParamVal();
    fitParams paramErr = GlobalEnv::instance()->DefaultParamErr();

    std::stringstream defaultparamsname;
    defaultparamsname << "defaultparams" << GlobalEnv::instance()->outputFileNameSuffix() << ".dat";
    std::ofstream theStreamDefault ( defaultparamsname.str().c_str() );

    GlobalEnv::instance()->fitParamsBase()->dumpParams(theStreamDefault, paramVal, paramErr);
}

void AppBase::generate(fitParams& theParams){
    std::shared_ptr<PwaGen> pwaGenPtr(new PwaGen());
    GlobalEnv::instance()->fitParamsBase()->printParams(theParams);
    pwaGenPtr->generate(GlobalEnv::instance()->Channel()->Lh(), theParams);
}

void AppBase::readEvents(EventList& theEventList, std::vector<std::string>& fileNames, ChannelID channelID, bool withEvtWeight, int evtStart, int evtStop){
  int noFinalStateParticles=GlobalEnv::instance()->Channel(channelID)->noFinalStateParticles();
  EventReaderDefault eventReader(fileNames, noFinalStateParticles, 0, withEvtWeight);
  eventReader.setUnit(GlobalEnv::instance()->parser()->unitInFile());
  eventReader.setOrder(GlobalEnv::instance()->parser()->orderInFile());

  if(GlobalEnv::instance()->Channel(channelID)->useMassRange()){
   eventReader.setMassRange(theEventList, GlobalEnv::instance()->Channel(channelID)->massRangeMin(),
                                          GlobalEnv::instance()->Channel(channelID)->massRangeMax(), GlobalEnv::instance()->Channel(channelID)->particleIndicesMassRange());
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

void AppBase::qaMode(fitParams& theStartParams, double evtWeightSumData, int noOfFreeFitParams){
  
  double theLh=GlobalEnv::instance()->Channel()->Lh()->calcLogLh(theStartParams);
  double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
  double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
  double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/(evtWeightSumData-noOfFreeFitParams-1);

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
      (new WaveContribution(GlobalEnv::instance()->Channel()->Lh(), theStartParams, thePwaCovMatrix));
  }
  else{
    theWaveContribution = std::shared_ptr<WaveContribution>
      (new WaveContribution(GlobalEnv::instance()->Channel()->Lh(), theStartParams));
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
}

void AppBase::qaModeSimple(EventList& dataEventList, EventList& mcEventList, fitParams& theStartParams, std::shared_ptr<EvtDataBaseList> evtDataBaseList, int noOfFreeFitParams){

  std::shared_ptr<AbsLh> absLh=GlobalEnv::instance()->Channel()->Lh();
  LHData theLHData;
  std::shared_ptr<WaveContribution> theWaveContribution(new WaveContribution(GlobalEnv::instance()->Channel()->Lh(), theStartParams));

  fitParams currentParams = theStartParams;

  for(int i=-1; i<static_cast<int>(theWaveContribution->NoOfContributions());i++){

    std::string contributionName="";

    if(i!=-1){
      contributionName = theWaveContribution->GetContributionName(i);
      MnUserParameters uPar = theWaveContribution->GetParametersForContribution(i);
      GlobalEnv::instance()->fitParamsBase()->getFitParamVal(uPar.Params(), currentParams);
    }
 
    std::shared_ptr<AbsHist> histPtr = GlobalEnv::instance()->Channel()->CreateHistInstance(contributionName);

    absLh->updateFitParams(currentParams);
  
    //loop over data events
    Event* anEvent;
    int evtCount = 0;
  
    dataEventList.rewind();
    while ((anEvent = dataEventList.nextEvent())){
      EvtData* currentDataEvt=evtDataBaseList->convertEvent(anEvent, evtCount);
      absLh->addDataToLogLh(currentDataEvt, currentParams, theLHData);
      histPtr->fillEvt(currentDataEvt, currentDataEvt->evtWeight, "data");
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

    if(i!=-1)
      continue;

    double theLh=absLh->mergeLogLhData(theLHData);
    double evtWeightSumData=theLHData.weightSum;
    double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
    double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
    double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/(evtWeightSumData-noOfFreeFitParams-1);
    double integralDataWoWeight=(double) dataEventList.size();
    
    std::ostringstream qaSummaryFileName;
    std::string outputFileNameSuffix= GlobalEnv::instance()->outputFileNameSuffix();
    qaSummaryFileName << "qaSummarySimple" << outputFileNameSuffix << ".dat";
    std::ofstream theQaStream ( qaSummaryFileName.str().c_str() );
    
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

    theQaStream.close();
  } // loop over contributions

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

void AppBase::fixParams(MnUserParameters& upar, std::vector<std::string> fixedParams){

  // Always fix the primary channel's scaling parameters
  std::string fixedScaleParam = GlobalEnv::instance()->Channel()->Lh()->getChannelScaleParam() + "Other";
  fixedParams.push_back(fixedScaleParam);
  Info << "Fixing scaling parameter " << fixedScaleParam << endmsg;

  const std::vector<MinuitParameter> theParams= upar.Parameters();
  std::vector<std::string> parNames;
  std::vector<MinuitParameter>::const_iterator itPar;
  for (itPar=theParams.begin(); itPar!=theParams.end(); ++itPar){
    parNames.push_back(itPar->GetName());
  }

  std::vector<std::string>::const_iterator itFix;
  for (itFix=fixedParams.begin(); itFix!=fixedParams.end(); ++itFix){
    //check if name exisists
    if(std::find(parNames.begin(), parNames.end(), (*itFix)) != parNames.end()) upar.Fix( (*itFix) );
    else{
      Alert << "parameter with name\t" << (*itFix) <<"\tdoes not exist!!!" << endmsg;
      exit(0);
    }
  }
}

void AppBase::fixAllReleaseScaleParams(MnUserParameters& upar){
  std::string scaleParam = GlobalEnv::instance()->Channel()->Lh()->getChannelScaleParam() + "Other";
  const std::vector<MinuitParameter> theParams= upar.Parameters();

  std::vector< MinuitParameter >::const_iterator itPar;
  for(itPar=theParams.begin(); itPar!=theParams.end(); ++itPar){
    //   if (scaleParam != itPar->GetName()) itPar->Fix();
    if (scaleParam != itPar->GetName()) upar.Fix( itPar->GetName() );
    else upar.Release( itPar->GetName() );
    // upar.Fix( itPar->GetName() );
    // if (scaleParam == itPar->GetName()) upar.Release( itPar->GetName() );
  }

}

FunctionMinimum AppBase::migradDefault(AbsFcn& theFcn, MnUserParameters& upar){
  MnMigrad migrad(theFcn, upar);
  Info <<"start migrad "<< endmsg;
  FunctionMinimum funcMin = migrad(0, GlobalEnv::instance()->parser()->tolerance());

  if(funcMin.IsValid()){
     return funcMin;
  }

  // Two more tries to get a valid result unsing strategy 2
  for(int j=0; j<2; j++){
     Warning <<"FM is invalid, try with strategy = 2."<< endmsg;

     // Check minimum covariance matrix
     bool badCovarianceDiagonal=false;
     if(funcMin.HasCovariance()){
	badCovarianceDiagonal = !PwaCovMatrix::DiagonalIsValid(funcMin.UserCovariance());
     }

     if(badCovarianceDiagonal){
       Warning << "Using default errors" << endmsg;
       MnUserParameters newParams = upar;
       for(unsigned int i=0; i< funcMin.UserParameters().Params().size();i++){
	  newParams.SetValue(i, funcMin.UserParameters().Params().at(i));
       }
       MnMigrad migrad2(theFcn, newParams, MnStrategy(2));
       funcMin = migrad2(0, GlobalEnv::instance()->parser()->tolerance());
    }
    else{
       MnUserParameters newParams = upar;
       for(unsigned int i=0; i< funcMin.UserParameters().Params().size();i++){
	  newParams.SetValue(i, funcMin.UserParameters().Params().at(i));
	  newParams.SetError(i, funcMin.UserParameters().Errors().at(i));
       }
       MnMigrad migrad2(theFcn, newParams, MnStrategy(2));
       funcMin = migrad2(0, GlobalEnv::instance()->parser()->tolerance());
    }

    if(funcMin.IsValid()){
       break;
    }
  }

  return funcMin;
}

void AppBase::printFitResult(FunctionMinimum& min, fitParams& theStartparams, std::ostream& os, double evtWeightSumData, int noOfFreeFitParams){

    double theLh = min.Fval();

    os << "\n\n********************** Final fit parameters *************************\n";
    os << "\n" << min.UserParameters() << "\n";
    os << "\n\n**************** Minuit FunctionMinimum information ******************" << std::endl;
    if(min.IsValid())             os << "\n Function minimum is valid.\n";
    else                          os << "\n WARNING: Function minimum is invalid!" << std::endl;
    if(min.HasValidCovariance())  os << "\n Covariance matrix is valid." << std::endl;
    else                          os << "\n WARNING: Covariance matrix is invalid!" << std::endl;
    os <<"\n Final LH: "<< std::setprecision(10) << theLh << "\n" << std::endl;
    os <<" # of function calls: " << min.NFcn() << std::endl;
    os <<" minimum edm: " << std::setprecision(10) << min.Edm()<<std::endl;
    if(!min.HasValidParameters()) os << " hasValidParameters() returned FALSE" << std::endl;
    if(!min.HasAccurateCovar())   os << " hasAccurateCovar() returned FALSE" << std::endl;
    if(!min.HasPosDefCovar()){    os << " hasPosDefCovar() returned FALSE" << std::endl;
                                  if(min.HasMadePosDefCovar()) os << " hasMadePosDefCovar() returned TRUE" << std::endl;
    }
    if(!min.HasCovariance())      os << " hasCovariance() returned FALSE" << std::endl;
    if(min.HasReachedCallLimit()) os << " hasReachedCallLimit() returned TRUE" << std::endl;
    if(min.IsAboveMaxEdm())       os << " isAboveMaxEdm() returned TRUE" << std::endl;
    if(min.HesseFailed())         os << " hesseFailed() returned TRUE" << std::endl;
    os << std::endl;

    MnUserParameters finalUsrParameters=min.UserParameters();
    const std::vector<double> finalParamVec=finalUsrParameters.Params();
    fitParams finalFitParams=theStartparams;
    GlobalEnv::instance()->fitParamsBase()->getFitParamVal(finalParamVec, finalFitParams);

    const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
    fitParams finalFitErrs=theStartparams;
    GlobalEnv::instance()->fitParamsBase()->getFitParamVal(finalParamErrorVec, finalFitErrs);

    std::ostringstream finalResultname;
    finalResultname << "finalResult" << GlobalEnv::instance()->outputFileNameSuffix() << ".dat";

    std::ofstream theStream ( finalResultname.str().c_str() );
    GlobalEnv::instance()->fitParamsBase()->dumpParams(theStream, finalFitParams, finalFitErrs);

    MnUserCovariance theCovMatrix = min.UserCovariance();
    std::ostringstream serializationFileName;
    serializationFileName << "serializedOutput" << GlobalEnv::instance()->outputFileNameSuffix() << ".dat";
    std::ofstream serializationStream(serializationFileName.str().c_str());
    boost::archive::text_oarchive boostOutputArchive(serializationStream);

    if(min.HasValidCovariance()){
        const PwaCovMatrix thePwaCovMatrix(theCovMatrix, finalUsrParameters, finalFitParams);
        boostOutputArchive << thePwaCovMatrix;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // calculate AIC, BIC criteria and output selected wave contrib
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Info <<"theLh = "<< theLh << endmsg;

    double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
    double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
    double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/(evtWeightSumData-noOfFreeFitParams-1);
    Info << "noOfFreeFitParams:\t" <<noOfFreeFitParams;
    Info << "evtWeightSumData:\t" <<evtWeightSumData;
    Info << "BIC:\t" << BICcriterion << endmsg;
    Info << "AIC:\t" << AICcriterion << endmsg;
    Info << "AICc:\t" << AICccriterion << endmsg;
}

bool AppBase::calcAndSendClientLh(NetworkClient& theClient, fitParams& theStartparams, ChannelID channelID){

  while(true){
    if(!theClient.WaitForParams()) return false;

    const std::vector<double> currentParamVec=theClient.GetParams();
    fitParams currentFitParams=theStartparams;
    GlobalEnv::instance()->fitParamsBase()->getFitParamVal(currentParamVec, currentFitParams);

    LHData theLHData;
    GlobalEnv::instance()->Channel(channelID)->Lh()->calcLogLhDataClient(currentFitParams, theLHData);

    if(!theClient.SendLH(theLHData.logLH_data, theLHData.LH_mc)) return false;
  }
  return true;
}

