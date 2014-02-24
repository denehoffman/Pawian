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

// HeliDecAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/HeliDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

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
    complex<double> amp = currentJPClamlam->parityFactor*theMag*expi*conj(theData->WignerDStringId.at(_wignerDKey).at(IdJLamXLam12));
    //    complex<double> amp = currentJPClamlam->parityFactor*theMag*expi*conj(theData->WignerDsString.at(_wignerDKey).at(currentJPClamlam->J).at(lamX).at(lambda));
    result+=amp;
  }
  //  result*=sqrt((2.*_JPCPtr->J+1.)/12.56637);
  result*=_preFactor*_isospinCG*sqrt(2.*_JPCPtr->J+1.);
  return result;
}




complex<double> HeliDecAmps::XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);

  if( fabs(lamX) > _JPCPtr->J) return result;

  int evtNo=theData->evtNo;
  Id2StringType currentSpinIndex=FunctionUtils::spin2Index(lamX,lamFs);

  if ( _cacheAmps && !_recalculate){
    result=_cachedAmpShortMap.at(evtNo).at(currentSpinIndex);
    //    result=_cachedAmpMap.at(evtNo).at(lamX).at(lamFs);
    result*=_absDyn->eval(theData, grandmaAmp);
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

    complex<double> amp = it->first->parityFactor*theMag*expi*conj( theData->WignerDStringId.at(_wignerDKey).at(IdJLamXLam12));
    result+=amp*daughterAmp(lambda1, lambda2, theData, lamFs);
  }

  result*=_preFactor*_isospinCG*sqrt(2.*_JPCPtr->J+1.);

  if ( _cacheAmps){
     theMutex.lock();
     _cachedAmpShortMap[evtNo][currentSpinIndex]=result;
     //     _cachedAmpMap[evtNo][lamX][lamFs]=result;
     theMutex.unlock();
}

  result*=_absDyn->eval(theData, grandmaAmp);
  if(result.real()!=result.real()){
    Alert << "result:\t" << result << endmsg;
    exit(0);
  }
  return result;
}

void  HeliDecAmps::getDefaultParams(fitParams& fitVal, fitParams& fitErr){

  std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
  std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
  std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
  std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

  std::vector< std::shared_ptr<const JPClamlam> >::const_iterator itlamlam;
  for(itlamlam=_JPClamlams.begin(); itlamlam!=_JPClamlams.end(); ++itlamlam){
    currentMagValMap[*itlamlam]=_factorMag;
    currentPhiValMap[*itlamlam]=0.;
    currentMagErrMap[*itlamlam]=_factorMag;
    currentPhiErrMap[*itlamlam]=0.3;
  }

  fitVal.MagLamLams[_key]=currentMagValMap;
  fitVal.PhiLamLams[_key]=currentPhiValMap;
  fitErr.MagLamLams[_key]=currentMagErrMap;
  fitErr.PhiLamLams[_key]=currentPhiErrMap;

  _absDyn->getDefaultParams(fitVal, fitErr);

  if(!_daughter1IsStable) _decAmpDaughter1->getDefaultParams(fitVal, fitErr);
  if(!_daughter2IsStable) _decAmpDaughter2->getDefaultParams(fitVal, fitErr);
}

void HeliDecAmps::print(std::ostream& os) const{
  return; //dummy
}


bool HeliDecAmps::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;

   if(_absDyn->checkRecalculation(theParamVal)) _recalculate=true;

   if(!_daughter1IsStable) {
     if(_decAmpDaughter1->checkRecalculation(theParamVal)) _recalculate=true;
   }
   if(!_daughter2IsStable){
     if(_decAmpDaughter2->checkRecalculation(theParamVal)) _recalculate=true;
   }

   if(!_recalculate){
     std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >& magMap=theParamVal.MagLamLams[_key];
     std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >& phiMap=theParamVal.PhiLamLams[_key];
     std::vector< std::shared_ptr<const JPClamlam> >::iterator it;
     for (it=_JPClamlams.begin(); it!=_JPClamlams.end(); ++it){
       double theMag=magMap[*it];
       double thePhi=phiMap[*it];

       if (!CheckDoubleEquality(theMag, _currentParamMagLamLams[*it])){
	 _recalculate=true;
	 return _recalculate;
       }
       if (!CheckDoubleEquality(thePhi, _currentParamPhiLamLams[*it])){
	 _recalculate=true;
	 return _recalculate;
       }
     }
   }
   return _recalculate;
}


void  HeliDecAmps::updateFitParams(fitParams& theParamVal){
   std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >& magMap=theParamVal.MagLamLams[_key];
   std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >& phiMap=theParamVal.PhiLamLams[_key];

   std::vector< std::shared_ptr<const JPClamlam> >::iterator it;
   for (it=_JPClamlams.begin(); it!=_JPClamlams.end(); ++it){
     double theMag=magMap[*it];
     double thePhi=phiMap[*it];
     _currentParamMagLamLams[*it]=theMag;
     _currentParamPhiLamLams[*it]=thePhi;

     std::vector< std::shared_ptr<const JPClamlam> >& currentLPClamlamVec=_JPClamlamSymMap[*it];
     std::vector< std::shared_ptr<const JPClamlam> >::iterator itLamLam;
     for (itLamLam=currentLPClamlamVec.begin(); itLamLam!=currentLPClamlamVec.end(); ++itLamLam){
       _currentParamMagLamLams[*itLamLam]=theMag;
       _currentParamPhiLamLams[*itLamLam]=thePhi;
     }
   }


   _absDyn->updateFitParams(theParamVal);

   if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(theParamVal);
   if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(theParamVal);

}
