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

// HeliDecAmps class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/HeliDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"
#include "FitParams/AbsPawianParameters.hh"

HeliDecAmps::HeliDecAmps(std::shared_ptr<IsobarHeliDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_JPClamlams(theDec->JPClamlamAmps())
  ,_factorMag(1.)
{
  if(_JPClamlams.size()>0) _factorMag=1./sqrt(_JPClamlams.size());
  Particle* daughter1=_decay->daughter1Part();
  Particle* daughter2=_decay->daughter2Part();
  _parityFactor=_JPCPtr->P*daughter1->theParity()*daughter2->theParity()*pow(-1,_JPCPtr->J-daughter1->J()-daughter2->J());
  Info << "_parityFactor=\t" << _parityFactor << endmsg;

  bool identicalDaughters=false;
  if( (*daughter1)==(*daughter2)) identicalDaughters=true;

  //fill JPClamlamSymMap
  std::vector< std::shared_ptr<const JPClamlam> >::iterator it;
  for(it=_JPClamlams.begin(); it!=_JPClamlams.end(); ++it){
    std::shared_ptr<const JPClamlam> currentSym(new JPClamlam(*it, -(*it)->lam1, -(*it)->lam2, _parityFactor));
    //    std::vector< std::shared_ptr<const JPClamlam> > currentLPClamlamVec=_JPClamlamSymMap[*it];
    _JPClamlamSymMap[*it].push_back(currentSym);
    if(identicalDaughters && ( fabs((*it)->lam1) != fabs((*it)->lam2) ) ){
    	std::shared_ptr<const JPClamlam> currentSymIdPart1(new JPClamlam(*it, (*it)->lam2, (*it)->lam1, (*it)->parityFactor));
    	std::shared_ptr<const JPClamlam> currentSymIdPart2(new JPClamlam(*it, -(*it)->lam2, -(*it)->lam1, _parityFactor));
    	_JPClamlamSymMap[*it].push_back(currentSymIdPart1);
    	_JPClamlamSymMap[*it].push_back(currentSymIdPart2);
    }
  }
}

HeliDecAmps::HeliDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
{
  Particle* daughter1=_decay->daughter1Part();
  Particle* daughter2=_decay->daughter2Part();
  _parityFactor=_JPCPtr->P*daughter1->theParity()*daughter2->theParity()*pow(-1,_JPCPtr->J-daughter1->J()-daughter2->J());
  Info << "_parityFactor=\t" << _parityFactor << endmsg;
}

HeliDecAmps::~HeliDecAmps()
{
}

void  HeliDecAmps::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){

  std::vector< std::shared_ptr<const JPClamlam> >::const_iterator itlamlam;
  for(itlamlam=_JPClamlams.begin(); itlamlam!=_JPClamlams.end(); ++itlamlam){
    //fill magnitude
    std::string magName=(*itlamlam)->name()+_key+"Mag";
    double valMag=_factorMag;
    double errMag=_factorMag/2.;
    // double minMag=0.;
    // double maxMag=_factorMag+30.*errMag;

    fitPar->Add(magName, valMag, errMag);
    //    fitPar->SetLimits(magName, minMag, maxMag);

    std::string phiName=(*itlamlam)->name()+_key+"Phi";
    double valPhi=0.;
    double errPhi=0.2;
    //no limits for phi parameter
    fitPar->Add(phiName, valPhi, errPhi);
  }

  _absDyn->fillDefaultParams(fitPar);

  if(!_daughter1IsStable) _decAmpDaughter1->fillDefaultParams(fitPar);
  if(!_daughter2IsStable) _decAmpDaughter2->fillDefaultParams(fitPar);
}

void HeliDecAmps::fillParamNameList(){
  _paramNameList.clear();
  std::vector< std::shared_ptr<const JPClamlam> >::const_iterator itlamlam;
  for(itlamlam=_JPClamlams.begin(); itlamlam!=_JPClamlams.end(); ++itlamlam){
    //magnitude
    std::string magName=(*itlamlam)->name()+_key+"Mag";
    _paramNameList.push_back(magName);

    //    Info << "HeliDecAmps add parameter " << magName << " to paramNameList" << endmsg;

    //phase
    std::string phiName=(*itlamlam)->name()+_key+"Phi";
    _paramNameList.push_back(phiName);

    //    Info << "HeliDecAmps add parameter " << phiName << " to paramNameList" << endmsg;
  }

}


void HeliDecAmps::print(std::ostream& os) const{
  return; //dummy
}

void HeliDecAmps::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){

  std::vector< std::shared_ptr<const JPClamlam> >::const_iterator itlamlam;
  for(itlamlam=_JPClamlams.begin(); itlamlam!=_JPClamlams.end(); ++itlamlam){

    std::string magName=(*itlamlam)->name()+_key+"Mag";
    std::string phiName=(*itlamlam)->name()+_key+"Phi";

    double theMag=fabs(fitPar->Value(magName));
    double thePhi=fitPar->Value(phiName);

    _currentParamMagLamLams[*itlamlam]=theMag;
    _currentParamPhiLamLams[*itlamlam]=thePhi;
    complex<double> expi(cos(thePhi), sin(thePhi));
    
    _currentParamMagExpi[*itlamlam]=theMag*expi;
    _currentParamPreFacMagExpi[*itlamlam]=_preFactor*_isospinCG*sqrt(2.*_JPCPtr->J+1.)*theMag*expi;

    std::vector< std::shared_ptr<const JPClamlam> >& currentLPClamlamSymVec=_JPClamlamSymMap.at(*itlamlam);
     std::vector< std::shared_ptr<const JPClamlam> >::iterator itLamLamSym;
     for (itLamLamSym=currentLPClamlamSymVec.begin(); itLamLamSym!=currentLPClamlamSymVec.end(); ++itLamLamSym){
       _currentParamMagLamLams[*itLamLamSym]=theMag;
       _currentParamPhiLamLams[*itLamLamSym]=thePhi;
       _currentParamMagExpi[*itLamLamSym]=theMag*expi;
       _currentParamPreFacMagExpi[*itLamLamSym]=_preFactor*_isospinCG*sqrt(2.*_JPCPtr->J+1.)*theMag*expi;
     }
  }

   _absDyn->updateFitParams(fitPar);

   if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(fitPar);
   if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(fitPar);
}

void HeliDecAmps::printCurrentAmpParams(Spin& lamX, Spin& lamFs){
  Info << "Amp name: " << name() << endmsg;
  std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator it;
    for(it=_currentParamMagLamLams.begin(); it!=_currentParamMagLamLams.end(); ++it){
      Spin lambda1= it->first->lam1;
      Spin lambda2= it->first->lam2;
      Spin lambda = lambda1-lambda2;
      if( fabs(lambda) > it->first->J) continue;
      
      if(_enabledlamFsDaughter1 && lamFs!=lambda1) continue;
      if(_enabledlamFsDaughter2 && lamFs!=lambda2) continue;
      it->first->print(std::cout);
      Info << " it->first->parityFactor: " << it->first->parityFactor << endmsg;
      Info << "_currentParamPreFacMagExpi.at(it->first): " << _currentParamPreFacMagExpi.at(it->first) << endmsg;
    }
}


complex<double> HeliDecAmps::XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){
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
    complex<double> amp = currentJPClamlam->parityFactor*theMag*expi*conj(theData->WignerDStringStringId.at(_wignerDKey).at(_decay->wignerDRefKey()).at(IdJLamXLam12));
    result+=amp;
  }

  result*=_preFactor*_isospinCG*sqrt(2.*_JPCPtr->J+1.);
  return result;
}


complex<double> HeliDecAmps::XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

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
    //    complex<double> amp = it->first->parityFactor*_currentParamPreFacMagExpi.at(it->first)*conj( theData->WignerDStringId.at(_wignerDKey).at(IdJLamXLam12));
    complex<double> amp = it->first->parityFactor*_currentParamPreFacMagExpi.at(it->first)*conj( theData->WignerDStringStringId.at(_wignerDKey).at(_decay->wignerDRefKey()).at(IdJLamXLam12) );

    //    complex<double> amp = it->first->parityFactor*_currentParamPreFacMagExpi.at(it->first)*conj( theData->WignerDStringStringId.at(_wignerDKey).at(_decay->wignerDRefKey()).at(IdJLamXLam12));
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

void HeliDecAmps::calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp){
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
