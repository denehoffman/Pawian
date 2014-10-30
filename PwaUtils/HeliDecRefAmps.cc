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

// HeliDecRefAmps class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/HeliDecRefAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

HeliDecRefAmps::HeliDecRefAmps(std::shared_ptr<IsobarHeliDecay> theDec, ChannelID channelID) :
  AbsHeliDecAmps(theDec, channelID)
{
}

HeliDecRefAmps::HeliDecRefAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsHeliDecAmps(theDec, channelID)
{
}

HeliDecRefAmps::~HeliDecRefAmps()
{
}


complex<double> HeliDecRefAmps::XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){
  complex<double> result(0.,0.);

  std::string refKey=_refKey;
  if (0!=grandmaAmp) refKey=grandmaAmp->refKey();


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
    complex<double> amp = currentJPClamlam->parityFactor*theMag*expi*conj(theData->WignerDStringStringId.at(_wignerDKey).at(refKey).at(IdJLamXLam12));
    result+=amp;
  }

  result*=_preFactor*_isospinCG*sqrt(2.*_JPCPtr->J+1.);
  return result;
}




complex<double> HeliDecRefAmps::XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);

  std::string refKey=_refKey;
  if (0!=grandmaAmp) refKey=grandmaAmp->refKey();

  if( fabs(lamX) > _JPCPtr->J) return result;

  int evtNo=theData->evtNo;
  Id2StringType currentSpinIndex=FunctionUtils::spin2Index(lamX,lamFs);

  if ( _cacheAmps && !_recalculate){
    result=_cachedAmpMapNew.at(evtNo).at(refKey).at(_absDyn->grandMaKey(grandmaAmp)).at(currentSpinIndex);
    //    result*=_absDyn->eval(theData, grandmaAmp);
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

    double theMag=it->second;
    double thePhi=_currentParamPhiLamLams.at(it->first);
    complex<double> expi(cos(thePhi), sin(thePhi));
    unsigned int IdJLamXLam12=FunctionUtils::spin3Index(_J, lamX, lambda);

    complex<double> amp = it->first->parityFactor*theMag*expi*conj( theData->WignerDStringStringId.at(_wignerDKey).at(refKey).at(IdJLamXLam12));
    result+=amp*daughterAmp(lambda1, lambda2, theData, lamFs);
  }

  result*=_preFactor*_isospinCG*sqrt(2.*_JPCPtr->J+1.);

  // if(absDec()->useProdBarrier()){
  //   result *= BarrierFactor::BlattWeisskopf(absDec()->orbMomMin(), theData->DoubleString.at(_wignerDKey), BarrierFactor::qRDefault) /
  //     BarrierFactor::BlattWeisskopf(absDec()->orbMomMin(), theData->DoubleString.at(_wignerDKey + "qNorm"), BarrierFactor::qRDefault);
  // }
  // else result*=_absDyn->eval(theData, grandmaAmp, absDec()->orbMomMin());

  result*=_absDyn->eval(theData, grandmaAmp, absDec()->orbMomMin());
  
  if(result.real()!=result.real()){
    Alert << "result:\t" << result << endmsg;
    exit(0);
  }

  if ( _cacheAmps){
    theMutex.lock();
    _cachedAmpMapNew[evtNo][refKey][_absDyn->grandMaKey(grandmaAmp)][currentSpinIndex]=result;
    theMutex.unlock();
  }

  return result;
}

