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

// HeliMultipoleDecNonRefAmps class definition file. -*- C++ -*-
// Copyright 2015 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <sstream>
#include <mutex>

#include "PwaUtils/HeliMultipoleDecNonRefAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"
#include "FitParams/AbsPawianParameters.hh"

HeliMultipoleDecNonRefAmps::HeliMultipoleDecNonRefAmps(std::shared_ptr<IsobarHeliDecay> theDec, ChannelID channelID) :
  HeliDecAmps(theDec, channelID)
  ,_noOfAmps(0)
{
  _daughterGamma=theDec->daughter1Part();
  if(_daughterGamma->twoJ() != 2 || _daughterGamma->theParity() != -1 ||  !_daughterGamma->theCParity() || _daughterGamma->mass() > 1.e-6){
    Alert << "daughter particle 1 is not a photon" << endmsg;
    Alert << "twoJ(): " << _daughterGamma->twoJ() << "\ttheParity(): " << _daughterGamma->theParity() << "\tCParity(): " << _daughterGamma->theCParity() <<"\tmass(): " << _daughterGamma->mass() << endmsg; 
    exit(0); 
  }

  _daughter2=theDec->daughter2Part();
  Spin JgammaMin(_JPCPtr->J-_daughter2->J());
  if(JgammaMin<0) JgammaMin=Spin(_daughter2->J()-_JPCPtr->J);
  if(JgammaMin==0) JgammaMin=Spin(1);

  int parityPiPf=_JPCPtr->P*_daughter2->theParity();

  Spin JgammaMax=Spin(_JPCPtr->J+_daughter2->J());
  
  _noOfAmps=_JPClamlams.size();
  _currentParamLocalMagExpi.resize(_noOfAmps);
  _MagParamNames.resize(_noOfAmps);
  _PhiParamNames.resize(_noOfAmps);
  _JgammaVec.resize(_noOfAmps);
  
  for(int i=0; i< _noOfAmps; ++i){
    Spin currentSpin(JgammaMin+i);
    _JgammaMap[i]=currentSpin;
    _JgammaVec.push_back(currentSpin);
    std::string currentKey=_key;
    std::stringstream tmpSpinStr;
    tmpSpinStr << currentSpin;
    int currentSpinInt=currentSpin;
     if(parityPiPf<0){
      if(int(currentSpinInt)%2 ==0 ) currentKey+="M"+tmpSpinStr.str();
      else currentKey+="E"+tmpSpinStr.str();
    }
    else{ //parityPiPf==1
      if(int(currentSpinInt)%2 ==0 ) currentKey+="E"+tmpSpinStr.str();
      else currentKey+="M"+tmpSpinStr.str();
    } 
    _MagParamNames[i]=currentKey+"_Mag";
    _PhiParamNames[i]=currentKey+"_Phi";
  }
}

HeliMultipoleDecNonRefAmps::HeliMultipoleDecNonRefAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  HeliDecAmps(theDec, channelID)
  ,_noOfAmps(0)
{
  _daughterGamma=theDec->daughter1Part();
  if(_daughterGamma->twoJ() != 2 || _daughterGamma->theParity() != -1 ||  !_daughterGamma->theCParity() || _daughterGamma->mass() > 1.e-6){
    Alert << "daughter particle 1 is not a photon" << endmsg;
    exit(0); 
  }

  _daughter2=theDec->daughter2Part();
  Spin JgammaMin(_JPCPtr->J-_daughter2->J());
  if(JgammaMin<0) JgammaMin=Spin(_daughter2->J()-_JPCPtr->J);
  if(JgammaMin==0) JgammaMin=Spin(1);

  int parityPiPf=_JPCPtr->P*_daughter2->theParity();

  Spin JgammaMax=Spin(_JPCPtr->J+_daughter2->J());
  
  _noOfAmps=_JPClamlams.size();
  _currentParamLocalMagExpi.resize(_noOfAmps);
  _MagParamNames.resize(_noOfAmps);
  _PhiParamNames.resize(_noOfAmps);
  _JgammaVec.resize(_noOfAmps);
  
  for(int i=0; i< _noOfAmps; ++i){
    Spin currentSpin(JgammaMin+i);
    _JgammaMap[i]=currentSpin;
    _JgammaVec.push_back(currentSpin);
    std::string currentKey=_key;
    std::stringstream tmpSpinStr;
    tmpSpinStr << currentSpin;
    int currentSpinInt=currentSpin;
    if(parityPiPf<0){
      if(currentSpinInt%2 ==0 ) currentKey+="M"+tmpSpinStr.str();
      else currentKey+="E"+tmpSpinStr.str();
    }
    else{ //parityPiPf==1
      if(currentSpinInt%2 ==0 ) currentKey+="E"+tmpSpinStr.str();
      else currentKey+="M"+tmpSpinStr.str();
    } 
    _MagParamNames[i]=currentKey+"_Mag";
    _PhiParamNames[i]=currentKey+"_Phi";
  }
}

HeliMultipoleDecNonRefAmps::~HeliMultipoleDecNonRefAmps()
{
}


void  HeliMultipoleDecNonRefAmps::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){

  for (int i=0; i<_noOfAmps; ++i){
    fitPar->Add(_MagParamNames.at(i), 1., 0.5);
    //    fitPar->SetLimits(_MagParamNames.at(i), 0., 20.);

    fitPar->Add(_PhiParamNames.at(i), 0., 0.2); 
  }

  _absDyn->fillDefaultParams(fitPar);
  if(!_daughter2IsStable) _decAmpDaughter2->fillDefaultParams(fitPar);
}

void HeliMultipoleDecNonRefAmps::fillParamNameList(){
  _paramNameList.clear();

  for (int i=0; i<_noOfAmps; ++i){
    _paramNameList.push_back(_MagParamNames.at(i));
    //    InfoMsg << "HeliMultipoleDecNonRefAmps: add parameter " << _MagParamNames.at(i) << " to paramNameList" << endmsg;
    _paramNameList.push_back(_PhiParamNames.at(i));
    //    InfoMsg << "HeliMultipoleDecNonRefAmps: add parameter " << _PhiParamNames.at(i) << " to paramNameList" << endmsg;
  }

}

void HeliMultipoleDecNonRefAmps::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  
  for (int i=0; i<_noOfAmps; ++i){
    double theLocalMag=fabs(fitPar->Value(_MagParamNames.at(i)));
    double theLocalPhi=fitPar->Value(_PhiParamNames.at(i));
    _currentParamLocalMagExpi[i]=std::polar(theLocalMag, theLocalPhi);
  }

  _absDyn->updateFitParams(fitPar);
  if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(fitPar);
}


complex<double> HeliMultipoleDecNonRefAmps::XdecAmp(Spin& lamX, EvtData* theData, AbsXdecAmp* grandmaAmp){
  complex<double> result(0.,0.);
  if( fabs(lamX) > _JPCPtr->J) return result;
 

  short currentSpinIndex=FunctionUtils::spin1IdIndex(_projIdThreadMap.at(std::this_thread::get_id()), lamX); 

  if (!_recalculate){
      result=_cachedAmpIdMap.at(theData->evtNo).at(_absDyn->grandMaId(grandmaAmp)).at(currentSpinIndex);
      if(result.real()!=result.real()) {
        DebugMsg << "result:\t" << result << endmsg;
      }
      return result;
    }

  //  std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator it;
  //  InfoMsg << "_currentParamMagLamLams.size(): " << _currentParamMagLamLams.size() << endmsg;
  
  //  for(it=_currentParamMagLamLams.begin(); it!=_currentParamMagLamLams.end(); ++it){

  std::vector< std::shared_ptr<const JPClamlam> >::iterator it;
  for(it=_JPClamlams.begin(); it!=_JPClamlams.end(); ++it){
    //Spin lambda1= it->first->lam1;  //gamma
    //Spin lambda2= it->first->lam2;  //X
    Spin lambda1= (*it)->lam1;  //gamma                                                         
    Spin lambda2= (*it)->lam2;  //X
    //    InfoMsg << "main: J: " << (*it)->J << " lamX: " << lamX
    //	    << " lam1: " << lambda1
    //	    << " lam2: " << lambda2 << endmsg;
    result+=heliAmpLoop(theData, lamX, lambda1, lambda2, (*it)->J, false);
    if(lambda1!=0 || lambda2!=0){
      lambda1=-lambda1;
      lambda2=-lambda2;
      result+=heliAmpLoop(theData, lamX, lambda1, lambda2, (*it)->J, true);
    }
  }

  if (_absDyn->isLdependent()) result*=_cachedDynLMap.at(std::this_thread::get_id());
   else result*=_cachedDynIdMap.at(std::this_thread::get_id()).at(_absDyn->grandMaId(grandmaAmp));

  if(result.real()!=result.real()){
    Alert << "result:\t" << result << endmsg;
    exit(0);
  }

  if ( _cacheAmps){
     _cachedAmpIdMap[theData->evtNo][_absDyn->grandMaId(grandmaAmp)][currentSpinIndex]=result;
  }

  return result;
}

complex<double> HeliMultipoleDecNonRefAmps::heliAmpLoop(EvtData* theData, Spin& lamX, Spin& lam1, Spin& lam2, const Spin& J, bool isSym){
  complex<double> result(0.,0.);
  Spin lambda = lam2-lam1;
  bool doCalc=true;
  if( fabs(lambda) > J) doCalc=false;
  // InfoMsg << "_daughter1IsStable: " << _daughter1IsStable << " _lam1MinThreadMap.at(std::this_thread::get_id()): " << _lam1MinThreadMap.at(std::this_thread::get_id()) << endmsg;
  
  if(_daughter1IsStable && (_lam1MinThreadMap.at(std::this_thread::get_id())!=lam1)) doCalc=false;
  if(_daughter2IsStable && (_lam2MinThreadMap.at(std::this_thread::get_id())!=lam2)) doCalc=false;
  //  InfoMsg << "J: " << J << " lamX: " << lamX << " lam1: " << lam1 << " lam2: " << lam2 << " lambda: " << lambda << " doCalc: " << doCalc << endmsg;
  if(doCalc){
   unsigned int IdJLamXLam12=FunctionUtils::spin3Index(_J, lamX, lambda);
   complex<double> currentAmp(0.,0.);
   for (int i=0; i<_noOfAmps; ++i){
      double currentParity=_daughter2->theParity();
      if(int(_JgammaMap.at(i))%2 == 0) currentParity *= -1.;
      double parityFactor=1.;
      if(currentParity<0. && isSym) parityFactor=-1.;
      //if(int(_JgammaMap.at(i))%2 == 0 && lam1<0) parityFactor *= -1.; //this is a fix which is needed to be checked
      currentAmp+=sqrt(2.*_JgammaMap.at(i)+1.)
        *parityFactor
        *Clebsch(_JgammaMap.at(i), -lam1, _daughter2->J(), lam2, _JPCPtr->J, lambda)
        *_currentParamLocalMagExpi.at(i);
  }

    result+=currentAmp*conj( theData->WignerDIdId3.at(_decay->wigDWigDRefId()).at(IdJLamXLam12) )*daughterAmp(lam1, lam2, theData);
  }
  return result;
}
