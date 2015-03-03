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

// LSDecNonRefAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/LSDecNonRefAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

LSDecNonRefAmps::LSDecNonRefAmps(std::shared_ptr<IsobarLSDecay> theDec, ChannelID channelID) :
  AbsLSDecAmps(theDec, channelID)
  ,_Smax(0)
{
  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=_LSs.begin(); it!=_LSs.end(); ++it){
    if( (*it)->S > _Smax ) _Smax=(*it)->S;
  }
}

LSDecNonRefAmps::LSDecNonRefAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsLSDecAmps(theDec, channelID)
  ,_Smax(0)
{
  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=_LSs.begin(); it!=_LSs.end(); ++it){
    if( (*it)->S > _Smax ) _Smax=(*it)->S;
  }
}

LSDecNonRefAmps::~LSDecNonRefAmps()
{
}


complex<double> LSDecNonRefAmps::XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);

  if(fixDaughterNr == 1){
     _lam1Min = _lam1Max = lamDec;
  }
  else if(fixDaughterNr == 2){
     _lam2Min = _lam2Max = lamDec;
  }
  else{
     Alert << "Invalid fixDaughterNr in XdecPartAmp." << endmsg;
  }

  if(_enabledlamFsDaughter1){
    _lam1Min=lamFs;
    _lam1Max=lamFs;
  }
  else if(_enabledlamFsDaughter2){
    _lam2Min=lamFs;
    _lam2Max=lamFs;
  }

  result=lsLoop( grandmaAmp, lamX, theData, _lam1Min, _lam1Max, _lam2Min, _lam2Max, false);
  //  result*=_preFactor*_isospinCG;
  //  if (!_absDyn->isLdependent()) result *=_absDyn->eval(theData, grandmaAmp);

  return result;
}




complex<double> LSDecNonRefAmps::XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);
  if( fabs(lamX) > _JPCPtr->J) return result;

  int evtNo=theData->evtNo;

  Id2StringType currentSpinIndex=FunctionUtils::spin2Index(lamX,lamFs); 
  //  unsigned short currentSpinIndex=lamX.ToIndex()*100+lamFs.ToIndex();

  if ( _cacheAmps && !_recalculate){
    result=_cachedAmpMap.at(evtNo).at(_absDyn->grandMaKey(grandmaAmp)).at(currentSpinIndex);
    return result;
  }

  if(_enabledlamFsDaughter1){
    _lam1Min=lamFs;
    _lam1Max=lamFs;
  }
  else if(_enabledlamFsDaughter2){
    _lam2Min=lamFs;
    _lam2Max=lamFs;
  }


  result=lsLoop(grandmaAmp, lamX, theData, _lam1Min, _lam1Max, _lam2Min, _lam2Max, true, lamFs);
  //  result*=_preFactor*_isospinCG;
  //  if (!_absDyn->isLdependent()) result *=_absDyn->eval(theData, grandmaAmp);

  if ( _cacheAmps){
     theMutex.lock();
     _cachedAmpMap[evtNo][_absDyn->grandMaKey(grandmaAmp)][currentSpinIndex]=result;
     theMutex.unlock();
  }

  if(result.real()!=result.real()){
    Info << "dyn name: " << _absDyn->name() 
	 << "\nname(): " << name()
	 << endmsg;
    Alert << "result:\t" << result << endmsg;
    exit(0);
  }
  return result;
}


complex<double> LSDecNonRefAmps::lsLoop(AbsXdecAmp* grandmaAmp, Spin& lamX, EvtData* theData, Spin& lam1Min, Spin& lam1Max, Spin& lam2Min, Spin& lam2Max, bool withDecs, Spin lamFs ){
 
  complex<double> result(0.,0.);

  std::vector< std::shared_ptr<const LScomb> >::iterator it;

  std::map<Id3StringType, complex<double> >& currentWignerDMap=theData->WignerDStringId.at(_wignerDKey);

  for(Spin lambda1=lam1Min; lambda1<=lam1Max; ++lambda1){

    for(Spin lambda2=lam2Min; lambda2<=lam2Max; ++lambda2){
      Spin lambda = lambda1-lambda2;
      if( fabs(lambda)>_JPCPtr->J || fabs(lambda)>_Smax) continue;

      map<std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& cgPre_LSMap=_cgPreFactor_LamLamLSMap.at(lambda1).at(lambda2);    
      complex<double> amp(0.,0.);     
      for (it=_LSs.begin(); it!=_LSs.end(); ++it){
	if( fabs(lambda)>(*it)->S) continue;
	if (_absDyn->isLdependent()) amp+=_currentParamPreFacMagExpi.at(*it)*cgPre_LSMap.at(*it)*_cachedDynLSMap.at(std::this_thread::get_id()).at((*it)->L);
	else amp+=_currentParamPreFacMagExpi.at(*it)*cgPre_LSMap.at(*it);
      }
      Id3StringType IdJLamXLam12=FunctionUtils::spin3Index(_J, lamX, lambda);
      amp *= conj(currentWignerDMap.at(IdJLamXLam12));
      if(withDecs) amp *=daughterAmp(lambda1, lambda2, theData, lamFs);
      result+=amp;    
    }
  }
  if (!_absDyn->isLdependent()) result *=_cachedDynMap.at(std::this_thread::get_id()).at(_absDyn->grandMaKey(grandmaAmp));
  return result;
}

void LSDecNonRefAmps::calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp){
  if(!_recalculate) return; 

  if(!_absDyn->isLdependent()){
    AbsXdecAmp::calcDynamics(theData, grandmaAmp);
    return;
  }

 std::vector< std::shared_ptr<const LScomb> >::iterator it;
 for (it=_LSs.begin(); it!=_LSs.end(); ++it){
   theMutex.lock();
   _cachedDynLSMap[std::this_thread::get_id()][(*it)->L]=_absDyn->eval(theData, grandmaAmp, (*it)->L);
   theMutex.unlock();
 }  

 if(!_daughter1IsStable) _decAmpDaughter1->calcDynamics(theData, this);
 if(!_daughter2IsStable) _decAmpDaughter2->calcDynamics(theData, this);
 return;
}

// void LSDecNonRefAmps::resetRecalcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp){
//   if(!_absDyn->isLdependent()){
//     AbsXdecAmp::resetRecalcDynamics(theData, grandmaAmp);
//     return;
//   }
//   _cachedDynRecalcLSMap[std::this_thread::get_id()][(*it)->L]=true;
//   if(!_daughter1IsStable) _decAmpDaughter1->resetRecalcDynamics(theData, this);
//   if(!_daughter2IsStable) _decAmpDaughter2->resetRecalcDynamics(theData, this);
//   return;
// }

