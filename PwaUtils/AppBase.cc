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
#include <fstream>
#include <string>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "PwaUtils/AppBase.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/AbsEnv.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/PwaGen.hh"
#include "PwaUtils/ParserBase.hh"
#include "PwaUtils/AbsHist.hh"
#include "PwaUtils/WaveContribution.hh"
#include "PwaUtils/PwaCovMatrix.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Event/Event.hh"
#include "Event/EventReaderDefault.hh"

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

