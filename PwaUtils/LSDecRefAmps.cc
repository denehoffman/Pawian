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

// LSDecRefAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/LSDecRefAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

LSDecRefAmps::LSDecRefAmps(std::shared_ptr<IsobarLSDecay> theDec, ChannelID channelID) :
  AbsLSDecAmps(theDec, channelID)
{
}

LSDecRefAmps::LSDecRefAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsLSDecAmps(theDec, channelID)
{
}

LSDecRefAmps::~LSDecRefAmps()
{
}


complex<double> LSDecRefAmps::XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){
  std::string refKey=_refKey;
  if (0!=grandmaAmp) refKey=grandmaAmp->refKey();

  complex<double> result(0.,0.);

  Spin lam1Min=-_Jdaughter1;
  Spin lam1Max= _Jdaughter1;
  Spin lam2Min=-_Jdaughter2;
  Spin lam2Max=_Jdaughter2;

  if(fixDaughterNr == 1){
     lam1Min = lam1Max = lamDec;
  }
  else if(fixDaughterNr == 2){
     lam2Min = lam2Max = lamDec;
  }
  else{
     Alert << "Invalid fixDaughterNr in XdecPartAmp." << endmsg;
  }

  if(_enabledlamFsDaughter1){
    lam1Min=lamFs;
    lam1Max=lamFs;
  }
  else if(_enabledlamFsDaughter2){
    lam2Min=lamFs;
    lam2Max=lamFs;
  }

  result=lsLoopRef(grandmaAmp, refKey, lamX, theData, lam1Min, lam1Max, lam2Min, lam2Max, false);

  return result;
}




complex<double> LSDecRefAmps::XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  std::string& refKey=_refKey;
  if (0!=grandmaAmp) refKey=grandmaAmp->refKey();
  
  complex<double> result(0.,0.);
  if( fabs(lamX) > _JPCPtr->J) return result;

  int evtNo=theData->evtNo;

  Id2StringType currentSpinIndex=FunctionUtils::spin2Index(lamX,lamFs); 

  if ( _cacheAmps && !_recalculate){
    result=_cachedAmpMapNew.at(evtNo).at(refKey).at(_absDyn->grandMaKey(grandmaAmp)).at(currentSpinIndex);
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


  result=lsLoopRef(grandmaAmp, refKey, lamX, theData, _lam1Min, _lam1Max, _lam2Min, _lam2Max, true, lamFs);

  if ( _cacheAmps){
     theMutex.lock();
     _cachedAmpMapNew[evtNo][refKey][_absDyn->grandMaKey(grandmaAmp)][currentSpinIndex]=result;
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


complex<double> LSDecRefAmps::lsLoopRef(AbsXdecAmp* grandmaAmp, std::string& refKey, Spin& lamX, EvtData* theData, Spin& lam1Min, Spin& lam1Max, Spin& lam2Min, Spin& lam2Max, bool withDecs, Spin lamFs ){
 
  complex<double> result(0.,0.);

  std::map<Id3StringType, complex<double> >& currentWignerDMap=theData->WignerDStringStringId.at(_wignerDKey).at(refKey);

  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=_LSs.begin(); it!=_LSs.end(); ++it){

    map<Spin,map<Spin, double > >& currentCgFactor=_cgPreFactor.at(*it);

    double theMag=_currentParamMags.at(*it);
    double thePhi=_currentParamPhis.at(*it);
    complex<double> expi(cos(thePhi), sin(thePhi));

    complex<double> tmpResult(0.,0.);
    for(Spin lambda1=lam1Min; lambda1<=lam1Max; ++lambda1){
      for(Spin lambda2=lam2Min; lambda2<=lam2Max; ++lambda2){
	Spin lambda = lambda1-lambda2;
	if( fabs(lambda)>_JPCPtr->J || fabs(lambda)>(*it)->S) continue;
	Id3StringType IdJLamXLam12=FunctionUtils::spin3Index(_J, lamX, lambda);
	complex<double> amp = theMag*expi*currentCgFactor.at(lambda1).at(lambda2)*conj(currentWignerDMap.at(IdJLamXLam12));

      	if(withDecs) amp *=daughterAmp(lambda1, lambda2, theData, lamFs);
	tmpResult+=amp;
      }
    }

    tmpResult*=_absDyn->eval(theData, grandmaAmp, (*it)->L);

    result+=tmpResult; 
  }

  result*=_preFactor*_isospinCG;
  // if(result.real()!=result.real()){
  //   Alert << "result:\t" << result << endmsg;
  //   exit(0);
  // }
  return result;
}


