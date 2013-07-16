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
#include "PwaUtils/AbsEnv.hh"
#include "PwaUtils/FitParamsBase.hh"
//#include "PwaUtils/AbsFcn.hh"
#include "PwaUtils/PwaGen.hh"
#include "PwaUtils/ParserBase.hh"
#include "PwaUtils/AbsHist.hh"
#include "PwaUtils/WaveContribution.hh"
#include "PwaUtils/PwaCovMatrix.hh"
#include "PwaUtils/NetworkClient.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Event/Event.hh"
#include "Event/EventReaderDefault.hh"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserCovariance.h"

AppBase::AppBase(AbsEnv* absEnv, std::shared_ptr<AbsLh> theLhPtr, std::shared_ptr<FitParamsBase> theFitParamBase) :
  _absEnv(absEnv), 
  _absLhPtr(theLhPtr),
  _fitParamBasePtr(theFitParamBase)
{
}

AppBase::~AppBase()
{
}

void AppBase::dumpDefaultParams(){
    fitParams defaultVal;
    fitParams defaultErr;
    _absLhPtr->getDefaultParams(defaultVal, defaultErr);

    std::stringstream defaultparamsname;
    defaultparamsname << "defaultparams" << _absEnv->outputFileNameSuffix() << ".dat";
    std::ofstream theStreamDefault ( defaultparamsname.str().c_str() );
    
    _fitParamBasePtr->dumpParams(theStreamDefault, defaultVal, defaultErr);
}

void AppBase::generate(fitParams& theParams){
    std::shared_ptr<PwaGen> pwaGenPtr(new PwaGen(_absEnv));
    pwaGenPtr->generate(_absLhPtr, theParams);
    _fitParamBasePtr->printParams(theParams);
}

void AppBase::readEvents(EventList& theEventList, std::vector<std::string>& fileNames, bool withEvtWeight, int evtStart, int evtStop){
  int noFinalStateParticles=_absEnv->noFinalStateParticles();
  EventReaderDefault eventReader(fileNames, noFinalStateParticles, 0, withEvtWeight);
  eventReader.setUnit(_absEnv->parser()->unitInFile());
  eventReader.setOrder(_absEnv->parser()->orderInFile());

  if(_absEnv->useMassRange())  eventReader.setMassRange(theEventList, _absEnv->massRangeMin(), _absEnv->massRangeMax(), _absEnv->particleIndicesMassRange());  
  
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

void AppBase::qaMode(fitParams& theStartParams, double evtWeightSumData, int noOfFreeFitParams ){

  double theLh=_absLhPtr->calcLogLh(theStartParams);
  Info <<"theLh = "<< theLh << endmsg;

  double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
  double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
  double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/(evtWeightSumData-noOfFreeFitParams-1);

  std::shared_ptr<WaveContribution> theWaveContribution;
  if(_absEnv->parser()->calcContributionError()){
    std::string serializationFileName = _absEnv->serializationFileName();
    std::ifstream serializationStream(serializationFileName.c_str());

    if(!serializationStream.is_open()){
      Alert << "Could not open serialization file." << endmsg;
      exit(0);
    }
       boost::archive::text_iarchive boostInputArchive(serializationStream);

       std::shared_ptr<PwaCovMatrix> thePwaCovMatrix(new PwaCovMatrix);
       boostInputArchive >> *thePwaCovMatrix;
       theWaveContribution = std::shared_ptr<WaveContribution>
	  (new WaveContribution(_absLhPtr, theStartParams, thePwaCovMatrix));
  }
  else{
    theWaveContribution = std::shared_ptr<WaveContribution>
      (new WaveContribution(_absLhPtr, theStartParams));
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
    std::string outputFileNameSuffix= _absEnv->outputFileNameSuffix();
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
}

void AppBase::fixParams(MnUserParameters& upar, const std::vector<std::string>& fixedParams){
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

FunctionMinimum AppBase::migradDefault(AbsFcn& theFcn, MnUserParameters& upar){
  MnMigrad migrad(theFcn, upar);
  Info <<"start migrad "<< endmsg;
  FunctionMinimum funcMin = migrad();
  if(!funcMin.IsValid()) {
    //try with higher strategy
    Info <<"FM is invalid, try with strategy = 2."<< endmsg;
    MnMigrad migrad2(theFcn, funcMin.UserState(), MnStrategy(2));
    funcMin = migrad2();
  }

  return funcMin;
}

void AppBase::printFitResult(FunctionMinimum& min, fitParams& theStartparams, std::ostream& os, std::string outputFileNameSuffix){

    os << "\n\n********************** Final fit parameters *************************\n";
    os << "\n" << min.UserParameters() << "\n";
    os << "\n\n**************** Minuit FunctionMinimum information ******************" << std::endl;
    if(min.IsValid())             os << "\n Function minimum is valid.\n";
    else                          os << "\n WARNING: Function minimum is invalid!" << std::endl;
    if(min.HasValidCovariance())  os << "\n Covariance matrix is valid." << std::endl;
    else                          os << "\n WARNING: Covariance matrix is invalid!" << std::endl;
    os <<"\n Final LH: "<< std::setprecision(10) << min.Fval() << "\n" << std::endl;
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
    _fitParamBasePtr->getFitParamVal(finalParamVec, finalFitParams);

    const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
    fitParams finalFitErrs=theStartparams;
    _fitParamBasePtr->getFitParamVal(finalParamErrorVec, finalFitErrs);

    std::ostringstream finalResultname;
    finalResultname << "finalResult" << outputFileNameSuffix << ".dat";

    std::ofstream theStream ( finalResultname.str().c_str() );
    _fitParamBasePtr->dumpParams(theStream, finalFitParams, finalFitErrs);

    MnUserCovariance theCovMatrix = min.UserCovariance();
    std::ostringstream serializationFileName;
    serializationFileName << "serializedOutput" << outputFileNameSuffix << ".dat";
    std::ofstream serializationStream(serializationFileName.str().c_str());
    boost::archive::text_oarchive boostOutputArchive(serializationStream);

    if(min.HasValidCovariance()){
        PwaCovMatrix thePwaCovMatrix(theCovMatrix, finalUsrParameters, finalFitParams);
        boostOutputArchive << thePwaCovMatrix;
    }

}

bool AppBase::calcAndSendClientLh(NetworkClient& theClient, fitParams& theStartparams){

    if(!theClient.WaitForParams()) return false;
    
    const std::vector<double> currentParamVec=theClient.GetParams();
    fitParams currentFitParams=theStartparams;
    _fitParamBasePtr->getFitParamVal(currentParamVec, currentFitParams);
    
    LHData theLHData;
    _absLhPtr->calcLogLhDataClient(currentFitParams, theLHData);
    
    if(!theClient.SendLH(theLHData.logLH_data, theLHData.weightSum, theLHData.LH_mc)) return false;
    return true;
}   

