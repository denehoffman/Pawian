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

// TensorPsiToGamXDecAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/TensorPsiToGamXDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/IsobarTensorPsiToGamXDecay.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

TensorPsiToGamXDecAmps::TensorPsiToGamXDecAmps(std::shared_ptr<IsobarTensorPsiToGamXDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_noOfAmps(0)
{
  _noOfAmps=theDec->noOfAmplitudes();
  _ampLMap=theDec->ampLMap();

  _currentParamLocalMags.resize(_noOfAmps);
  _currentParamLocalPhis.resize(_noOfAmps);
  _MagParamNames.resize(_noOfAmps);
  _PhiParamNames.resize(_noOfAmps);
  _MagParamNames[0]=_key+"Mag1";
  _PhiParamNames[0]=_key+"_1Phi";
  if(_noOfAmps>1){
    _MagParamNames[1]=_key+"Mag2";
    _PhiParamNames[1]=_key+"_2Phi";
  }
  if(_noOfAmps>2){
    _MagParamNames[2]=_key+"Mag3";
    _PhiParamNames[2]=_key+"_3Phi";
  }

}

TensorPsiToGamXDecAmps::TensorPsiToGamXDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_noOfAmps(0)
{
  //_noOfAmps=???
}

TensorPsiToGamXDecAmps::~TensorPsiToGamXDecAmps()
{
}


complex<double> TensorPsiToGamXDecAmps::XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){
  
  complex<double> result(0.,0.);
  
  // Spin lam1Min=-_Jdaughter1;
  // Spin lam1Max= _Jdaughter1;
  // Spin lam2Min=-_Jdaughter2;
  // Spin lam2Max=_Jdaughter2;

  // if(fixDaughterNr == 1){
  //    lam1Min = lam1Max = lamDec;
  // }
  // else if(fixDaughterNr == 2){
  //    lam2Min = lam2Max = lamDec;
  // }
  // else{
  //    Alert << "Invalid fixDaughterNr in XdecPartAmp." << endmsg;
  // }

  // if(_enabledlamFsDaughter1){
  //   lam1Min=lamFs;
  //   lam1Max=lamFs;
  // }
  // else if(_enabledlamFsDaughter2){
  //   lam2Min=lamFs;
  //   lam2Max=lamFs;
  // }

  //  complex<double> result=lsLoop(grandmaAmp, lamX, theData, lam1Min, lam1Max, lam2Min, lam2Max, false);

  return result;
}




complex<double> TensorPsiToGamXDecAmps::XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){
 
  complex<double> result(0.,0.);
  if( fabs(lamX) > _JPCPtr->J) return result;

  int evtNo=theData->evtNo;
  Id2StringType currentSpinIndex=FunctionUtils::spin2Index(lamX,lamFs);
  
  if ( _cacheAmps && !_recalculate){
    result=_cachedAmpMap.at(evtNo).at(_absDyn->grandMaKey(grandmaAmp)).at(currentSpinIndex);
    //    result*=_absDyn->eval(theData, grandmaAmp);
    return result;
  }

  Spin lam2Min=-_Jdaughter2;
  Spin lam2Max=_Jdaughter2;
  
  for(Spin lambda2=lam2Min; lambda2<=lam2Max; ++lambda2){
    complex<double> tmpResult(0.,0.);
    for(int i=0; i<_noOfAmps; ++i){
      double theMag=_currentParamLocalMags.at(i);
      double thePhi=_currentParamLocalPhis.at(i);
      complex<double> expi(cos(thePhi), sin(thePhi));
      tmpResult+=theMag*expi*theData->ComplexDoubleInt3SpinString.at(_name).at(i).at(lamX).at(lamFs).at(lambda2)*_absDyn->eval(theData, grandmaAmp,_ampLMap.at(i));
    }

    result+=tmpResult*daughterAmp(lambda2, theData, lamFs);    
  }

  result*=_absDyn->eval(theData, grandmaAmp);

  if ( _cacheAmps){
     theMutex.lock();
     _cachedAmpMap[evtNo][_absDyn->grandMaKey(grandmaAmp)][currentSpinIndex]=result;
     theMutex.unlock();
  }

  return result;
}

void  TensorPsiToGamXDecAmps::getDefaultParams(fitParams& fitVal, fitParams& fitErr){

  for (int i=0; i<_noOfAmps; ++i){
    fitVal.otherParams[_MagParamNames.at(i)]=15.;
    fitErr.otherParams[_MagParamNames.at(i)]=5.;
    fitVal.otherParams[_PhiParamNames.at(i)]=0.;
    fitErr.otherParams[_PhiParamNames.at(i)]=0.3;
  }

  _absDyn->getDefaultParams(fitVal, fitErr);

  if(!_daughter2IsStable) _decAmpDaughter2->getDefaultParams(fitVal, fitErr);
}

void TensorPsiToGamXDecAmps::print(std::ostream& os) const{
  return; //dummy
}


bool TensorPsiToGamXDecAmps::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;

   if(_absDyn->checkRecalculation(theParamVal)) _recalculate=true;

   if(!_daughter2IsStable){
     if(_decAmpDaughter2->checkRecalculation(theParamVal)){
       _recalculate=true;
       return _recalculate;
     }
   }

   if(!_recalculate){
     for (int i=0; i<_noOfAmps; ++i){
       double theMag=theParamVal.otherParams[_MagParamNames.at(i)];
       if(!CheckDoubleEquality(theMag, _currentParamLocalMags.at(i))){
	 _recalculate=true;
	 return _recalculate;
       }

       double thePhi=theParamVal.otherParams[_PhiParamNames.at(i)];
       if(!CheckDoubleEquality(thePhi, _currentParamLocalPhis.at(i))){
	 _recalculate=true;
	 return _recalculate;
       }
     }
   }

   return _recalculate;
}


void  TensorPsiToGamXDecAmps::updateFitParams(fitParams& theParamVal){

  for (int i=0; i<_noOfAmps; ++i){
    _currentParamLocalMags[i]=theParamVal.otherParams.at(_MagParamNames.at(i));
    _currentParamLocalPhis[i]=theParamVal.otherParams.at(_PhiParamNames.at(i));
  }
  _absDyn->updateFitParams(theParamVal);
  
  if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(theParamVal);

}

complex<double> TensorPsiToGamXDecAmps::daughterAmp(Spin& lam2, EvtData* theData, Spin& lamFs){
  complex<double> result(1.,0.);
  if(!_daughter2IsStable) result *= _decAmpDaughter2->XdecAmp(lam2, theData, lamFs, this);
  return result;
}

