//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// pipiScatteringBaseLh class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "pipiScatteringUtils/pipiScatteringBaseLh.hh"
#include "pipiScatteringUtils/PiPiScatteringChannelEnv.hh"
#include "ConfigParser/ParserBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/EvtDataScatteringList.hh"
#include "PwaUtils/PiPiScatteringXdecAmps.hh"
#include "ErrLogger/ErrLogger.hh"

#include "Utils/IdStringMapRegistry.hh"
#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


pipiScatteringBaseLh::pipiScatteringBaseLh(ChannelID channelID) :
  AbsLh(channelID)
{
  initialize();
}


pipiScatteringBaseLh::~pipiScatteringBaseLh()
{}

double pipiScatteringBaseLh::calcLogLh(std::shared_ptr<AbsPawianParameters> fitPar){
  _calcCounter++;
  updateFitParams(fitPar);

  double chi2All=0.; 
  //workaround for Tcheck
  if(GlobalEnv::instance()->Channel(_channelID)->parser()->productionFormalism() == "Tcheck"){
    _PiPiScatteringXdecAmp->doTcheck();
    return chi2All;
  }

  //loop over all data points
  std::vector<EvtData*>::iterator it;
  for(it=_evtDataVec.begin(); it!=_evtDataVec.end(); ++it){
    chi2All+=calcEvtIntensity( (*it), fitPar);
  }
  return chi2All;
}

complex<double> pipiScatteringBaseLh::calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData){

  complex<double> result(0.,0.);
  return result;
}


complex<double> pipiScatteringBaseLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData,
					     std::map <std::shared_ptr<const JPCLS>,
					     std::vector< std::shared_ptr<AbsXdecAmp> >,
					     pawian::Collection::SharedPtrLess > resAmps){
   complex<double> resultAmp(0.,0.);
   return resultAmp;
}

double pipiScatteringBaseLh::calcEvtIntensity(EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar){
  // chi2 fit must be included here
  // each point represents one event???

  Spin dummylamX(0);
  _PiPiScatteringXdecAmp->XdecAmp(dummylamX, theData);

  double dataPoint=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::DATA_PIPISCAT_NAME));
  double dataPointError=theData->DoubleId.at(IdStringMapRegistry::instance()->
					     stringId(EvtDataScatteringList::DATAERR_PIPISCAT_NAME));
  double fitPoint=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME));
  fitPoint *= fitPar->Value(_channelScaleParam);
  double chi2 = (dataPoint-fitPoint)*(dataPoint-fitPoint)/(dataPointError*dataPointError);
  double result=0.5*chi2;
  return result; 
}

void pipiScatteringBaseLh::print(std::ostream& os) const {

}


void  pipiScatteringBaseLh::initialize() {
  AbsLh::initialize();
  const std::shared_ptr<PiPiScatteringChannelEnv> pipiScatteringEnv =
    std::static_pointer_cast<PiPiScatteringChannelEnv>(GlobalEnv::instance()->PiPiScatteringChannel(_channelID));

  std::vector<std::shared_ptr<AbsDecay> > absDecList=pipiScatteringEnv->absDecayList()->getList();
  _XdecAmp=XdecAmpRegistry::instance()->getXdecAmp(_channelID, absDecList.at(0)->absDecPtr());
  _decAmps.push_back(_XdecAmp);
  _PiPiScatteringXdecAmp=std::dynamic_pointer_cast<PiPiScatteringXdecAmps>(_XdecAmp);
}


void pipiScatteringBaseLh::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  _XdecAmp->fillDefaultParams(fitPar);
  fitPar->Add(_channelScaleParam, 1., 0.01);
  fitPar->SetLimits(_channelScaleParam, 0., 3.);
}

void pipiScatteringBaseLh::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  _XdecAmp->updateFitParams(fitPar);
}

double pipiScatteringBaseLh::addDataToLogLh(EvtData* dataEvt, std::shared_ptr<AbsPawianParameters> fitPar, 
					    LHData& theLHData) {
  double intensity=calcEvtIntensity(dataEvt, fitPar);
  theLHData.logLH_data+=intensity;
  theLHData.weightSum+= dataEvt->evtWeight;
  return intensity;
}

double pipiScatteringBaseLh::addMcToLogLh(EvtData* mcEvt, std::shared_ptr<AbsPawianParameters> fitPar, LHData& theLHData) {
  double intensity=0.;
  theLHData.LH_mc+=intensity;
  theLHData.num_mc++;
  return intensity;
}

void pipiScatteringBaseLh::calcLogLhDataClient(std::shared_ptr<AbsPawianParameters> fitPar, LHData& theLHData) {
  _calcCounter++;
  if (_cacheAmps && _calcCounter>1) checkRecalculation(fitPar, _oldFitPar);
  updateFitParams(fitPar);

  int numData = _evtDataVec.size();

  int eventStepData = numData / _noOfThreads;

  std::vector<std::thread> theThreads;
  std::vector<LHData> threadDataVec;
  threadDataVec.resize(_noOfThreads);

  for(int i = 0; i<_noOfThreads;i++){
     int eventMin = i*eventStepData;
     int eventMax = (i==_noOfThreads-1) ? (_evtDataVec.size() - 1) : (i+1)*eventStepData - 1;

     theThreads.push_back(std::thread(&pipiScatteringBaseLh::ThreadfuncData, this, eventMin, eventMax,
                                      fitPar, std::ref(threadDataVec.at(i))));
  }
  for(auto it = theThreads.begin(); it != theThreads.end(); ++it) {
     (*it).join();
  }

  theThreads.clear();
  for(auto it = threadDataVec.begin(); it!= threadDataVec.end(); ++it) {
     theLHData.logLH_data += (*it).logLH_data;
     theLHData.weightSum += (*it).weightSum;
     theLHData.LH_mc += 0.;
     theLHData.num_mc += 0.;
  }

   if(_calcCounter<2) _oldFitPar = std::shared_ptr<AbsPawianParameters>(fitPar->Clone());
   else _oldFitPar->SetAllValues(fitPar->Params());
}

void pipiScatteringBaseLh::ThreadfuncData(unsigned int minEvent, unsigned int maxEvent,
					  std::shared_ptr<AbsPawianParameters> fitPar, LHData& theLHData) {
  for (unsigned int i=minEvent; i<=maxEvent; ++i) {
    addDataToLogLh(_evtDataVec.at(i), fitPar, theLHData);
  }
}

void pipiScatteringBaseLh::ThreadfuncMc(unsigned int minEvent, unsigned int maxEvent,
					std::shared_ptr<AbsPawianParameters> fitPar, LHData& theLHData) {
  return;
}

double pipiScatteringBaseLh::calcFitVal( EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar) {
  _XdecAmp->updateFitParams(fitPar);
  Spin dummylamX(0);
  _PiPiScatteringXdecAmp->XdecAmp(dummylamX, theData);
  double fitPoint=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME));
  return fitPoint;
}


