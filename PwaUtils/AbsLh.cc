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

// AbsLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <vector>
#include <thread>

#include <boost/thread.hpp>

#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/FsParticleProjections.hh"
#include "GlobalEnv.hh"
#include "ConfigParser/ParserBase.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

AbsLh::AbsLh(std::shared_ptr<AbsLh> theAbsLhPtr):
  AbsParamHandler()
  , _channelID(theAbsLhPtr->getChannelID())
  ,_evtDataVec(theAbsLhPtr->getDataVec())
  ,_evtMCVec(theAbsLhPtr->getMcVec())
  ,_usePhasespace(GlobalEnv::instance()->Channel(_channelID)->parser()->usePhaseSpaceHyp())
  ,_useCohPhasespace(GlobalEnv::instance()->Channel(_channelID)->parser()->useCohPhaseSpaceHyp())
  ,_CohPhasespaceKey("CohPhasespace")
  ,_calcCounter(0)
  ,_noOfThreads(GlobalEnv::instance()->parser()->noOfThreads())
{
  if(_noOfThreads>1){
    WarningMsg << "multiple threading not supperted anymore; changed noOfTheads from: "
	       << _noOfThreads << " to: " << 1 << " !!!" << endmsg;
    _noOfThreads=1;
  }
   std::stringstream stringStrChannelId;
   stringStrChannelId << _channelID;
   _phasespaceKey="Phasespace"+stringStrChannelId.str();
   _CohPhasespaceKey="CohPhasespace"+stringStrChannelId.str();
   initialize();
}

AbsLh::AbsLh(ChannelID channelID) :
  AbsParamHandler()
  , _channelID(channelID)
  ,_usePhasespace(GlobalEnv::instance()->Channel(_channelID)->parser()->usePhaseSpaceHyp())
  ,_useCohPhasespace(GlobalEnv::instance()->Channel(_channelID)->parser()->useCohPhaseSpaceHyp())
  ,_calcCounter(0)
  ,_noOfThreads(GlobalEnv::instance()->parser()->noOfThreads())
{
    if(_noOfThreads>1){
    WarningMsg << "multiple threading not supperted anymore; changed noOfTheads from: "
               << _noOfThreads << " to: " << 1 << " !!!" << endmsg;
    _noOfThreads=1;
  }
  std::stringstream stringStrChannelId;
   stringStrChannelId << _channelID;
   _phasespaceKey="Phasespace"+stringStrChannelId.str();
   _CohPhasespaceKey="CohPhasespace"+stringStrChannelId.str();
   initialize();
}

AbsLh::~AbsLh()
{
}



void AbsLh::initialize(){
    std::vector<Particle*> fsParticles=GlobalEnv::instance()->Channel(_channelID)->finalStateParticles();
    std::vector<Particle*>::iterator itParticle;

    std::ostringstream channelScaleParamStream;
    channelScaleParamStream << "channelType" << GlobalEnv::instance()->Channel(_channelID)->channelType();
    channelScaleParamStream << GlobalEnv::instance()->Channel(_channelID)->channelSubTypeName();
    if(GlobalEnv::instance()->Channel(_channelID)->parser()->scalingWChannelID())    channelScaleParamStream << "Id" << _channelID;    
    for (itParticle=fsParticles.begin(); itParticle != fsParticles.end(); ++itParticle){
       channelScaleParamStream << (*itParticle)->name();
    }
    channelScaleParamStream << "Scaling"+GlobalEnv::instance()->Channel(_channelID)->parser()->channelScalingSuffix();
    if(GlobalEnv::instance()->Channel(_channelID)->parser()->addChannelScalingId()){
      channelScaleParamStream << _channelID;
    }
    _channelScaleParam = channelScaleParamStream.str();
    _fsParticleProjections=GlobalEnv::instance()->Channel(_channelID)->getFsParticleProjectionsPtr();
}



double AbsLh::calcLogLh(std::shared_ptr<AbsPawianParameters> fitPar){
  _calcCounter++;
  if (_cacheAmps && _calcCounter>1) checkRecalculation(fitPar, _oldFitPar);
  updateFitParams(fitPar);

  LHData theLHData;

  for (unsigned int i=0; i<_evtDataVec.size(); ++i){
    addDataToLogLh(_evtDataVec.at(i), fitPar, theLHData);
  }

  for (unsigned int i=0; i<_evtMCVec.size(); ++i){
    addMcToLogLh(_evtMCVec.at(i), fitPar, theLHData);
  }

  if(_calcCounter<2) _oldFitPar = std::shared_ptr<AbsPawianParameters>(fitPar->Clone());
  else _oldFitPar->SetAllValues(fitPar->Params());
  return mergeLogLhData(theLHData);
}

double AbsLh::addDataToLogLh(EvtData* dataEvt, std::shared_ptr<AbsPawianParameters> fitPar, LHData& theLHData){
  double intensity=calcEvtIntensity(dataEvt, fitPar);
  theLHData.logLH_data+=(dataEvt->evtWeight)*log(intensity);
  theLHData.weightSum+= dataEvt->evtWeight;
  theLHData.squaredWeightSum+= dataEvt->evtWeight*dataEvt->evtWeight;
  return intensity;
}

double AbsLh::addMcToLogLh(EvtData* mcEvt, std::shared_ptr<AbsPawianParameters> fitPar, LHData& theLHData){
  double intensity=calcEvtIntensity(mcEvt, fitPar);
  theLHData.LH_mc+=intensity;
  theLHData.num_mc++;
  return intensity;
}

void AbsLh::calcLogLhDataClient(std::shared_ptr<AbsPawianParameters> fitPar, LHData& theLHData){
  _calcCounter++;
  if (_cacheAmps && _calcCounter>1) checkRecalculation(fitPar, _oldFitPar);
  updateFitParams(fitPar);

  for (unsigned int i=0; i<_evtDataVec.size(); ++i){
    addDataToLogLh(_evtDataVec.at(i), fitPar, theLHData);
  }

  for (unsigned int i=0; i<_evtMCVec.size(); ++i){
    addMcToLogLh(_evtMCVec.at(i), fitPar, theLHData);
  }

   if(_calcCounter<2) _oldFitPar = std::shared_ptr<AbsPawianParameters>(fitPar->Clone());
   else _oldFitPar->SetAllValues(fitPar->Params());
}

double AbsLh::mergeLogLhData(LHData& theLHData, ChannelID channelId){

  double logLH=0.;
  if(GlobalEnv::instance()->Channel(channelId)->channelType() == AbsChannelEnv::CHANNEL_PIPISCATTERING){
    logLH=theLHData.logLH_data;
  }
  else{
    double logLH_mc_Norm=0.;
    
    if (theLHData.LH_mc>0.) logLH_mc_Norm=log(theLHData.LH_mc/theLHData.num_mc);
    logLH=0.5*theLHData.weightSum *(theLHData.LH_mc/theLHData.num_mc-1.)*(theLHData.LH_mc/theLHData.num_mc-1.)
      -theLHData.logLH_data
      +theLHData.weightSum*logLH_mc_Norm;
    logLH=(theLHData.weightSum/theLHData.squaredWeightSum)*logLH;
  }
  return GlobalEnv::instance()->Channel(channelId)->parser()->nllScalingFactor()*logLH;
}



void AbsLh::setHyps( const std::map<const std::string, bool>& theMap, bool& theHyp, std::string& theKey){

  std::map<const std::string, bool>::const_iterator iter= theMap.find(theKey);

  if (iter !=theMap.end()){
    theHyp= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << theHyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << theKey << " does not exist!!!" <<endmsg;
    exit(0);
  }
}

void AbsLh::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){

  if(_usePhasespace){
    fitPar->Add(_phasespaceKey, 0.01, 0.05);
    fitPar->SetLimits(_phasespaceKey, 0., 3.);     
  }

  if(_useCohPhasespace){
    fitPar->Add(_CohPhasespaceKey+"Mag", 0.01, 0.05);
    fitPar->SetLimits(_CohPhasespaceKey+"Mag", 0., 3.);

    fitPar->Add(_CohPhasespaceKey+"Phi", 0., 0.2);
  }

  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itDecs;
  for(itDecs=_decAmps.begin(); itDecs!=_decAmps.end(); ++itDecs){
    (*itDecs)->fillDefaultParams(fitPar);
  }

  fitPar->Add(_channelScaleParam, 1., 0.01);
  fitPar->SetLimits(_channelScaleParam, 0., 20.);
}

void AbsLh::cacheAmplitudes(){
  _cacheAmps=true;
  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator it;
  for (it=_decAmps.begin(); it!=_decAmps.end(); ++it){
    (*it)->cacheAmplitudes();
  }
}

void AbsLh::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
std::vector< std::shared_ptr<AbsXdecAmp> >::iterator it;
  for (it=_decAmps.begin(); it!=_decAmps.end(); ++it){
    (*it)->updateFitParams(fitPar);
  }
}

bool AbsLh::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){
  bool result=true;
  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator it;
  for (it=_decAmps.begin(); it!=_decAmps.end(); ++it){
    if(!(*it)->checkRecalculation(fitParNew, fitParOld)) result=false;
  }
  return result;
}

void AbsLh::setDataVec(std::vector<EvtData*> theVec) {
  if(_evtDataVec.size()>0){
    Alert << "data vector already set!!!" << endmsg;
    exit(0);
  }

  _evtDataVec=theVec;
}

void AbsLh::setMcVec(std::vector<EvtData*> theVec) {
  if(_evtMCVec.size()>0){
    Alert << "mc vector already set!!!" << endmsg;
    exit(0);
  }

  _evtMCVec=theVec;
}



