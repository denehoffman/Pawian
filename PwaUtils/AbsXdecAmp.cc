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

// AbsXdecAmp class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/DynRegistry.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/FsParticleProjections.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"


AbsXdecAmp::AbsXdecAmp(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsParamHandler()
  , _channelID(channelID)
  , _decay(theDec)
  , _name(theDec->name())
  ,_JPCPtr(theDec->motherJPC())
  //  , _jpcDecsName(_JPCPtr->jpcname()+"To"+theDec->daughter1Part()->name()+"_"+theDec->daughter2Part()->name())
  , _isoKey(_JPCPtr->jpcname()+theDec->fitParSuffix())
  ,_absDyn(theDec->getDynamics())
  ,_fsParticleProjections(GlobalEnv::instance()->Channel(channelID)->getFsParticleProjectionsPtr())
  ,_daughter1ProjId(10000)
  ,_daughter2ProjId(10000)
  ,_isWeakDecay(theDec->isWeakDeacy())
  ,_isospinCG(theDec->isospinCG())
  ,_preFactor(theDec->preFactor())
  ,_key("_"+theDec->fitParSuffix())
  ,_wignerDKey(theDec->wignerDKey())
  //  ,_refKey(theDec->refKey())
  ,_daughter1IsStable(theDec->isDaughter1Stable())
  ,_daughter2IsStable(theDec->isDaughter2Stable())
  ,_enabledlamFsDaughter1(false)
  ,_enabledlamFsDaughter2(false)
  ,_J(theDec->motherJPC()->J)
  ,_isFormationAmp(false)
{
  if(theDec->type() =="FormationDecay"){
    _isFormationAmp=true;
    _jpcDecsName=_JPCPtr->jpcname()+"To"+theDec->daughter1Part()->name();
  }
  else{
    _jpcDecsName=_JPCPtr->jpcname()+"To"+theDec->daughter1Part()->name()+"_"+theDec->daughter2Part()->name();
  }
}

AbsXdecAmp::~AbsXdecAmp()
{
}

void AbsXdecAmp::initialize(){
  if(!_decay->useIsospin()) _isospinCG=1.;

  _absDyn = DynRegistry::instance()->getDynamics(_decay);

  if(!_daughter1IsStable){
    std::shared_ptr<AbsDecay> decDaughter1=_decay->decDaughter1();
    std::shared_ptr<AbsDynamics> currentDyn=decDaughter1->getDynamics();
    if (0!=currentDyn) decDaughter1->getDynamics()->addGrandMa(_decay);
    else{
      WarningMsg <<"no dynamics for\t" << decDaughter1->name() << endmsg;
    }
    _decAmpDaughter1=XdecAmpRegistry::instance()->getXdecAmp(_channelID, decDaughter1);
  }

  if(!_daughter2IsStable){
    std::shared_ptr<AbsDecay> decDaughter2=_decay->decDaughter2();
    
    std::shared_ptr<AbsDynamics> currentDyn=decDaughter2->getDynamics();
    if (0!=currentDyn) decDaughter2->getDynamics()->addGrandMa(_decay);
    else{
      WarningMsg <<"no dynamics for\t" << decDaughter2->name() << endmsg;
    }
    _decAmpDaughter2=XdecAmpRegistry::instance()->getXdecAmp(_channelID, decDaughter2);
  }


  _Jdaughter1=(Spin) _decay->daughter1Part()->J();
  _Jdaughter2=(Spin) _decay->daughter2Part()->J();
   
  _lam1Min = -_Jdaughter1;
  _lam1Max = _Jdaughter1;
  _lam2Min = -_Jdaughter2;
  _lam2Max = _Jdaughter2;

  _daughter1Name=_decay->daughter1Part()->name();
  _daughter2Name=_decay->daughter2Part()->name();

  if( _daughter1IsStable) _daughter1ProjId=_fsParticleProjections->fsParticleId(_daughter1Name);
  if( _daughter1IsStable && _Jdaughter1>0) _enabledlamFsDaughter1=true;
  if( _daughter2IsStable) _daughter2ProjId=_fsParticleProjections->fsParticleId(_daughter2Name);
  if( _daughter2IsStable && _Jdaughter2>0) _enabledlamFsDaughter2=true;
 
}

complex<double> AbsXdecAmp::daughterAmp(Spin& lam1, Spin& lam2, EvtData* theData){
  complex<double> result(1.,0.);
  if(!_daughter1IsStable) result *= _decAmpDaughter1->XdecAmp(lam1, theData, this);
  if(!_daughter2IsStable) result *= _decAmpDaughter2->XdecAmp(lam2, theData, this);
  return result;
}

void AbsXdecAmp::cacheAmplitudes(){
  _cacheAmps=true;
  _absDyn->cacheAmplitudes();
  if(!_daughter1IsStable) _decAmpDaughter1->cacheAmplitudes();
  if(!_daughter2IsStable) _decAmpDaughter2->cacheAmplitudes();  
}

bool AbsXdecAmp::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){
  _recalculate=false;
  if(_absDyn->checkRecalculation(fitParNew, fitParOld)) _recalculate=true;
  
  if(!_daughter1IsStable) {
    if(_decAmpDaughter1->checkRecalculation(fitParNew, fitParOld)) _recalculate=true;
  }
  
  if(!_daughter2IsStable){
    if(_decAmpDaughter2->checkRecalculation(fitParNew, fitParOld)) _recalculate=true;
  }
  
  if (!_recalculate) _recalculate=AbsParamHandler::checkRecalculation(fitParNew, fitParOld);

  return _recalculate;
}

void AbsXdecAmp::calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp){
  if(!_recalculate) return;

  if(!_absDyn->isLdependent()){
    theMutex.lock();
    //    _cachedDynMap[std::this_thread::get_id()][_absDyn->grandMaKey(grandmaAmp)] = _absDyn->eval( theData, grandmaAmp);
    _cachedDynIdMap[std::this_thread::get_id()][_absDyn->grandMaId(grandmaAmp)] = _absDyn->eval( theData, grandmaAmp);
    theMutex.unlock();
  }
  if(!_daughter1IsStable) _decAmpDaughter1->calcDynamics(theData, this);

  if(!_isFormationAmp){
    if(!_daughter2IsStable) _decAmpDaughter2->calcDynamics(theData, this);
  }

  return;
}

void AbsXdecAmp::setSpinProjections(int projId){
  std::vector<Spin> projections=_fsParticleProjections->spinProjections().at(projId);
  theMutex.lock();
  _projIdThreadMap[std::this_thread::get_id()]=projId;
  
  _lam1MinThreadMap[std::this_thread::get_id()]=_lam1Min;
  _lam1MaxThreadMap[std::this_thread::get_id()]=_lam1Max;
  _lam2MinThreadMap[std::this_thread::get_id()]=_lam2Min;
  _lam2MaxThreadMap[std::this_thread::get_id()]=_lam2Max;
  theMutex.unlock();

  if(_daughter1IsStable){
    Spin currentProjection=projections.at(_daughter1ProjId);
    theMutex.lock();
    _lam1MinThreadMap[std::this_thread::get_id()]=currentProjection;
    _lam1MaxThreadMap[std::this_thread::get_id()]=currentProjection;
    DebugMsg << "set spin projection: " << _daughter1Name << " _lam1Min: " << _lam1MinThreadMap[std::this_thread::get_id()] << " _lam1Max: " << _lam1MaxThreadMap[std::this_thread::get_id()] << endmsg;
    theMutex.unlock();
  }
  else _decAmpDaughter1->setSpinProjections(projId);
  
  if(!_isFormationAmp){ 
    if(_daughter2IsStable){
      Spin currentProjection=projections.at(_daughter2ProjId);
      theMutex.lock();
      _lam2MinThreadMap[std::this_thread::get_id()]=currentProjection;
      _lam2MaxThreadMap[std::this_thread::get_id()]=currentProjection;
      DebugMsg << "set spin projection: " << _daughter2Name << " _lam1Min: " << _lam2MinThreadMap[std::this_thread::get_id()] << " _lam1Max: " << _lam2MaxThreadMap[std::this_thread::get_id()] << endmsg;
      theMutex.unlock();
    }
    else _decAmpDaughter2->setSpinProjections(projId);
  }

  return;
}


