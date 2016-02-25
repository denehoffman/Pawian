//************************************************************************//
//									  //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// HeliDecNonRefAmps class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/HeliDecNonRefAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

HeliDecNonRefAmps::HeliDecNonRefAmps(std::shared_ptr<IsobarHeliDecay> theDec, ChannelID channelID) :
  AbsHeliDecAmps(theDec, channelID)
{
}

HeliDecNonRefAmps::HeliDecNonRefAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsHeliDecAmps(theDec, channelID)
{
}

HeliDecNonRefAmps::~HeliDecNonRefAmps()
{
}


complex<double> HeliDecNonRefAmps::XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){
  complex<double> result(0.,0.);

  bool lamFs_daughter1=false;
  if( _daughter1IsStable && _Jdaughter1>0) lamFs_daughter1=true;

  bool lamFs_daughter2=false;
  if( _daughter2IsStable && _Jdaughter2>0) lamFs_daughter2=true;

  std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator it;

  for(it=_currentParamMagLamLams.begin(); it!=_currentParamMagLamLams.end(); ++it){
    std::shared_ptr<const JPClamlam> currentJPClamlam=it->first;
    if( fabs(lamX) > currentJPClamlam->J) continue;

    Spin lambda1= currentJPClamlam->lam1;
    Spin lambda2= currentJPClamlam->lam2;
    Spin lambda = lambda1-lambda2;
    if( fabs(lambda) > currentJPClamlam->J) continue;
    if(lamFs_daughter1 && lamFs!=lambda1) continue;
    if(lamFs_daughter2 && lamFs!=lambda2) continue;
    if(fixDaughterNr==1 && lamDec!=lambda1) continue;
    if(fixDaughterNr==2 && lamDec!=lambda2) continue;

    double theMag=it->second;
    double thePhi=_currentParamPhiLamLams[currentJPClamlam];
    complex<double> expi(cos(thePhi), sin(thePhi));
    Id3StringType IdJLamXLam12=FunctionUtils::spin3Index(_J, lamX, lambda);
    complex<double> amp = currentJPClamlam->parityFactor*theMag*expi*conj(theData->WignerDStringId.at(_wignerDKey).at(IdJLamXLam12));
    result+=amp;
  }

  result*=_preFactor*_isospinCG*sqrt(2.*_JPCPtr->J+1.);
  return result;
}




complex<double> HeliDecNonRefAmps::XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);

  if( fabs(lamX) > _JPCPtr->J) return result;

  int evtNo=theData->evtNo;
  Id2StringType currentSpinIndex=FunctionUtils::spin2Index(lamX,lamFs);

  if ( _cacheAmps && !_recalculate){
    result=_cachedAmpMap.at(evtNo).at(_absDyn->grandMaKey(grandmaAmp)).at(currentSpinIndex);
    if(result.real()!=result.real()) DebugMsg << "result:\t" << result << endmsg;
    return result;
  }

  std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator it;

  for(it=_currentParamMagLamLams.begin(); it!=_currentParamMagLamLams.end(); ++it){

    Spin lambda1= it->first->lam1;
    Spin lambda2= it->first->lam2;
    Spin lambda = lambda1-lambda2;
    if( fabs(lambda) > it->first->J) continue;

    if(_enabledlamFsDaughter1 && lamFs!=lambda1) continue;
    if(_enabledlamFsDaughter2 && lamFs!=lambda2) continue;

    unsigned int IdJLamXLam12=FunctionUtils::spin3Index(_J, lamX, lambda);
    complex<double> amp = it->first->parityFactor*_currentParamPreFacMagExpi.at(it->first)*conj( theData->WignerDStringId.at(_wignerDKey).at(IdJLamXLam12));
    result+=amp*daughterAmp(lambda1, lambda2, theData, lamFs);
  }

  //  result*=_preFactor*_isospinCG*sqrt(2.*_JPCPtr->J+1.);

  if (_absDyn->isLdependent()) result*=_cachedDynLMap.at(std::this_thread::get_id());
  else result*=_cachedDynMap.at(std::this_thread::get_id()).at(_absDyn->grandMaKey(grandmaAmp));
  
  if(result.real()!=result.real()){
    Alert << "result:\t" << result << endmsg;
    printCurrentAmpParams(lamX, lamFs); 
    exit(0);
  }

  if ( _cacheAmps){
    theMutex.lock();
    _cachedAmpMap[evtNo][_absDyn->grandMaKey(grandmaAmp)][currentSpinIndex]=result;
    theMutex.unlock();
  }

  return result;
}

void HeliDecNonRefAmps::calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp){
  if(!_recalculate) return; 

  if(!_absDyn->isLdependent()){
    AbsXdecAmp::calcDynamics(theData, grandmaAmp);
    return;
  }

  theMutex.lock();
  _cachedDynLMap[std::this_thread::get_id()]=_absDyn->eval(theData, grandmaAmp, absDec()->orbMomMin());
  theMutex.unlock();

 if(!_daughter1IsStable) _decAmpDaughter1->calcDynamics(theData, this);
 if(!_daughter2IsStable) _decAmpDaughter2->calcDynamics(theData, this);
 return;
}

