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

HeliMultipoleDecNonRefAmps::HeliMultipoleDecNonRefAmps(std::shared_ptr<IsobarHeliDecay> theDec, ChannelID channelID) :
  HeliDecNonRefAmps(theDec, channelID)
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

  int parityPiPf=_JPCPtr->P*_daughter2->theCParity();

  Spin JgammaMax=Spin(_JPCPtr->J+_daughter2->J());
  
  _noOfAmps=_JPClamlams.size();
  _currentParamLocalMags.resize(_noOfAmps);
  _currentParamLocalPhis.resize(_noOfAmps);
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
  HeliDecNonRefAmps(theDec, channelID)
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

  int parityPiPf=_JPCPtr->P*_daughter2->theCParity();

  Spin JgammaMax=Spin(_JPCPtr->J+_daughter2->J());
  
  _noOfAmps=_JPClamlams.size();
  _currentParamLocalMags.resize(_noOfAmps);
  _currentParamLocalPhis.resize(_noOfAmps);
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


void HeliMultipoleDecNonRefAmps::getDefaultParams(fitParCol& fitVal, fitParCol& fitErr){

  for (int i=0; i<_noOfAmps; ++i){
    fitVal.otherParams[_MagParamNames.at(i)]=1.;
    fitErr.otherParams[_MagParamNames.at(i)]=0.5;
    fitVal.otherParams[_PhiParamNames.at(i)]=0.;
    fitErr.otherParams[_PhiParamNames.at(i)]=0.3;
  }

  _absDyn->getDefaultParams(fitVal, fitErr);
  if(!_daughter2IsStable) _decAmpDaughter2->getDefaultParams(fitVal, fitErr);
}

bool HeliMultipoleDecNonRefAmps::checkRecalculation(fitParCol& theParamVal){

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


void  HeliMultipoleDecNonRefAmps::updateFitParams(fitParCol& theParamVal){

  for (int i=0; i<_noOfAmps; ++i){
    double theLocalMag=theParamVal.otherParams.at(_MagParamNames.at(i));
    _currentParamLocalMags[i]=theLocalMag;
    double theLocalPhi=theParamVal.otherParams.at(_PhiParamNames.at(i));
    _currentParamLocalPhis[i]=theLocalPhi;
    
    complex<double> expiLocal(cos(theLocalPhi), sin(theLocalPhi));
    _currentParamLocalMagExpi[i]=theLocalMag*expiLocal;
  }

  std::vector< std::shared_ptr<const JPClamlam> >::iterator it;
   for (it=_JPClamlams.begin(); it!=_JPClamlams.end(); ++it){
     Spin lamGamma=Spin(1);
     Spin lam2=(*it)->lam2;
     if (lam2<0) lam2=-(*it)->lam2;
     _currentParamMagExpi[*it]=complex<double> (0.,0.);
     for (int i=0; i<_noOfAmps; ++i){
       _currentParamMagExpi[*it]+= sqrt(2.*_JgammaMap.at(i)+1.)*Clebsch(_JgammaMap.at(i), lamGamma, _JPCPtr->J, lam2-1, _daughter2->J(), lam2)*_currentParamLocalMagExpi.at(i);
     }

     _currentParamPreFacMagExpi[*it]=_preFactor*_currentParamMagExpi[*it];
     _currentParamMagLamLams[*it]=1.; //dummy
     _currentParamPhiLamLams[*it]=0.; //dummy     
     std::vector< std::shared_ptr<const JPClamlam> >& currentLPClamlamVec=_JPClamlamSymMap.at(*it);
     std::vector< std::shared_ptr<const JPClamlam> >::iterator itLamLam;
     for (itLamLam=currentLPClamlamVec.begin(); itLamLam!=currentLPClamlamVec.end(); ++itLamLam){
       _currentParamMagExpi[*itLamLam]=_currentParamMagExpi.at(*it);
       _currentParamPreFacMagExpi[*itLamLam]=_preFactor*_currentParamMagExpi.at(*it);
       _currentParamMagLamLams[*itLamLam]=1.; //dummy
       _currentParamPhiLamLams[*itLamLam]=0.; //dummy
     }
   }


   _absDyn->updateFitParams(theParamVal);

   if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(theParamVal);
   if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(theParamVal);

}
