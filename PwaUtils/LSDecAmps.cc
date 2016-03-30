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

// LSDecAmps class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/LSDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"
#include "FitParams/AbsPawianParameters.hh"

LSDecAmps::LSDecAmps(std::shared_ptr<IsobarLSDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_LSs(theDec->LSAmps())
  ,_factorMag(1.)
  ,_lam1Min(-_Jdaughter1)
  ,_lam1Max(_Jdaughter1)
  ,_lam2Min(-_Jdaughter2)
  ,_lam2Max(_Jdaughter2)
  ,_Smax(0)
{
  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=_LSs.begin(); it!=_LSs.end(); ++it){
    if( (*it)->S > _Smax ) _Smax=(*it)->S;
  }

  if(_LSs.size()>0) _factorMag=1./sqrt(_LSs.size());
  Particle* daughter1=_decay->daughter1Part();
  Particle* daughter2=_decay->daughter2Part();
  _parityFactor=daughter1->theParity()*daughter2->theParity()*pow(-1,_JPCPtr->J-daughter1->J()-daughter2->J());
  Info << "_parityFactor=\t" << _parityFactor << endmsg;
  fillCgPreFactor();
}

LSDecAmps::LSDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_lam1Min(-_Jdaughter1)
  ,_lam1Max(_Jdaughter1)
  ,_lam2Min(-_Jdaughter2)
  ,_lam2Max(_Jdaughter2)
  ,_Smax(0)
{
  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=_LSs.begin(); it!=_LSs.end(); ++it){
    if( (*it)->S > _Smax ) _Smax=(*it)->S;
  }

  Particle* daughter1=_decay->daughter1Part();
  Particle* daughter2=_decay->daughter2Part();
  _parityFactor=daughter1->theParity()*daughter2->theParity()*pow(-1,_JPCPtr->J-daughter1->J()-daughter2->J());
  Info << "_parityFactor=\t" << _parityFactor << endmsg;
  fillCgPreFactor();
}

LSDecAmps::~LSDecAmps()
{
}

complex<double> LSDecAmps::XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

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

  return result;
}

complex<double> LSDecAmps::XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);
  if( fabs(lamX) > _JPCPtr->J) return result;

  int evtNo=theData->evtNo;

  Id2StringType currentSpinIndex=FunctionUtils::spin2Index(lamX,lamFs); 

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

complex<double> LSDecAmps::lsLoop(AbsXdecAmp* grandmaAmp, Spin& lamX, EvtData* theData, Spin& lam1Min, Spin& lam1Max, Spin& lam2Min, Spin& lam2Max, bool withDecs, Spin lamFs ){
 
  complex<double> result(0.,0.);

  std::vector< std::shared_ptr<const LScomb> >::iterator it;

  //std::map<Id3StringType, complex<double> >& currentWignerDMap=theData->WignerDStringId.at(_wignerDKey);
  std::map<Id3StringType, complex<double> >& currentWignerDMap=theData->WignerDStringStringId.at(_wignerDKey).at(_decay->wignerDRefKey());

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
      //      amp *= conj(_currentWignerDMap->at(IdJLamXLam12));
      //      amp *= conj(_currentWignerDMap.at(IdJLamXLam12));
      if(withDecs) amp *=daughterAmp(lambda1, lambda2, theData, lamFs);
      result+=amp;    
    }
  }
  if (!_absDyn->isLdependent()) result *=_cachedDynMap.at(std::this_thread::get_id()).at(_absDyn->grandMaKey(grandmaAmp));
  return result;
}

void LSDecAmps::calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp){
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

void  LSDecAmps::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){

  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    //fill magnitude
    std::string magName=(*itLS)->name()+_key+"Mag";
    double valMag=_factorMag;
    double errMag=_factorMag/2.;
    //    double minMag=0.;
    //    double maxMag=_factorMag+30.*errMag;

    fitPar->Add(magName, valMag, errMag);
    //    fitPar->SetLimits(magName, minMag, maxMag);

    std::string phiName=(*itLS)->name()+_key+"Phi";
    double valPhi=0.;
    double errPhi=0.2;
    //no limits for phi parameter
    fitPar->Add(phiName, valPhi, errPhi);
  }

  _absDyn->fillDefaultParams(fitPar);

  if(!_daughter1IsStable) _decAmpDaughter1->fillDefaultParams(fitPar);
  if(!_daughter2IsStable) _decAmpDaughter2->fillDefaultParams(fitPar);
}

void LSDecAmps::fillParamNameList(){
  _paramNameList.clear();
  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    //magnitude
    std::string magName=(*itLS)->name()+_key+"Mag";
    _paramNameList.push_back(magName);
    //phi
    std::string phiName=(*itLS)->name()+_key+"Phi";
    _paramNameList.push_back(phiName);
  }

}

void LSDecAmps::print(std::ostream& os) const{
  return; //dummy
}


void LSDecAmps::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){

  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    //fill magnitude
    std::string magName=(*itLS)->name()+_key+"Mag";
    double theMag= fabs(fitPar->Value(magName));
    
    std::string phiName=(*itLS)->name()+_key+"Phi";
    double thePhi=fitPar->Value(phiName);

    _currentParamMags[*itLS]=theMag;
    _currentParamPhis[*itLS]=thePhi;
    complex<double> expi(cos(thePhi), sin(thePhi));
    _currentParamMagExpi[*itLS]=theMag*expi;
    _currentParamPreFacMagExpi[*itLS]=_preFactor*_isospinCG*theMag*expi;
  }

   _absDyn->updateFitParams(fitPar);

  if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(fitPar);
  if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(fitPar);
}

void  LSDecAmps::fillCgPreFactor(){

  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=_LSs.begin(); it!=_LSs.end(); ++it){
    for(Spin lambda1=-_Jdaughter1; lambda1<=_Jdaughter1; ++lambda1){
      for(Spin lambda2=-_Jdaughter2; lambda2<=_Jdaughter2; ++lambda2){
	Spin lambda = lambda1-lambda2;
	if( fabs(lambda)>_JPCPtr->J || fabs(lambda)>(*it)->S) continue;

	_cgPreFactor[*it][lambda1][lambda2]=sqrt(2.*(*it)->L+1)
	  *Clebsch((*it)->L, 0, (*it)->S, lambda, _JPCPtr->J, lambda)
	  *Clebsch(_Jdaughter1, lambda1, _Jdaughter2, -lambda2, (*it)->S, lambda  );

	_cgPreFactor_LamLamLSMap[lambda1][lambda2][*it]=sqrt(2.*(*it)->L+1)
	  *Clebsch((*it)->L, 0, (*it)->S, lambda, _JPCPtr->J, lambda)
	  *Clebsch(_Jdaughter1, lambda1, _Jdaughter2, -lambda2, (*it)->S, lambda  );
      }
    }
  }
}

// void LSDecAmps::retrieveWignerDs(EvtData* theData){
//   //  _currentWignerDMap = &theData->WignerDStringStringId.at(_wignerDKey).at(_decay->wignerDRefKey());
//   // _currentWignerDMap = theData->WignerDStringStringId.at(_wignerDKey).at(_decay->wignerDRefKey());
// }
