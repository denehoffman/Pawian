//************************************************************************//
//									  //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// HeliDecAmps class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

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
  ,_isWeakDecay(false)
{
  initialize();
  if(_JPClamlams.size()>0) _factorMag=1./sqrt(_JPClamlams.size());
  Particle* daughter1=_decay->daughter1Part();
  Particle* daughter2=_decay->daughter2Part();

  //look whether it is a weak decay or not
  Particle* motherParticle=_decay->motherPart();
  if(0!=motherParticle){
    if( motherParticle->strange() != (daughter1->strange()+daughter2->strange())) _isWeakDecay=true;
    else if ( motherParticle->charm() != (daughter1->charm()+daughter2->charm())) _isWeakDecay=true;     
  }

  if(_isWeakDecay){
    _isospinCG=1.;
    _parityFactor=1.;
    InfoMsg << "decay " << _name << " is a weak decay; parity conservation etc. disabled!!!!" << endmsg;

    unsigned int JPClamlamsSize=_JPClamlams.size();
    for (unsigned int i=0; i<JPClamlamsSize; ++i){
      std::shared_ptr<const JPClamlam> currentJPClamlam=_JPClamlams.at(i);
      if(currentJPClamlam->lam1 == 0 && currentJPClamlam->lam2 == 0) continue; 
      std::shared_ptr<const JPClamlam> toBeAddedJPClamlam(new JPClamlam(currentJPClamlam, 
									-currentJPClamlam->lam1, 
									-currentJPClamlam->lam2, 
									_parityFactor));
      _JPClamlams.push_back(toBeAddedJPClamlam);      
    }
  }

  else{
    _parityFactor=_JPCPtr->P*daughter1->theParity()*daughter2->theParity()
      * pow(-1,_JPCPtr->J-daughter1->J()-daughter2->J());
    InfoMsg << "_parityFactor=\t" << _parityFactor << endmsg;
    
    bool identicalDaughters=false;
    if( (*daughter1)==(*daughter2)) identicalDaughters=true;
    
    //fill JPClamlamSymMap
    std::vector< std::shared_ptr<const JPClamlam> >::iterator it;
    for(it=_JPClamlams.begin(); it!=_JPClamlams.end(); ++it){
      if( (*it)->lam1 == 0 && (*it)->lam2 == 0) continue;
      std::shared_ptr<const JPClamlam> currentSym(new JPClamlam(*it, -(*it)->lam1, -(*it)->lam2, _parityFactor));
      _JPClamlamSymMap[*it].push_back(currentSym);
      if(identicalDaughters && ( fabs((*it)->lam1) != fabs((*it)->lam2) ) ){
	std::shared_ptr<const JPClamlam> currentSymIdPart1(new JPClamlam(*it, (*it)->lam2, (*it)->lam1, 1));
    	std::shared_ptr<const JPClamlam> currentSymIdPart2(new JPClamlam(*it, -(*it)->lam2, -(*it)->lam1, _parityFactor));
    	_JPClamlamSymMap[*it].push_back(currentSymIdPart1);
    	_JPClamlamSymMap[*it].push_back(currentSymIdPart2);
      }
    }
  }
}

HeliDecAmps::HeliDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_isWeakDecay(false)
{
  initialize();
  Particle* daughter1=_decay->daughter1Part();
  Particle* daughter2=_decay->daughter2Part();
  _parityFactor=_JPCPtr->P*daughter1->theParity()*daughter2->theParity()
    * pow(-1,_JPCPtr->J-daughter1->J()-daughter2->J());
  InfoMsg << "_parityFactor=\t" << _parityFactor << endmsg;
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
    double errMag=_factorMag/4.;

    fitPar->Add(magName, valMag, errMag);

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

    //phase
    std::string phiName=(*itlamlam)->name()+_key+"Phi";
    _paramNameList.push_back(phiName);
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

    complex<double> heliAmpCompl=std::polar(fabs(fitPar->Value(magName)), fitPar->Value(phiName));
    complex<double> currentHeliAmp=_preFactor*_isospinCG*sqrt(2.*_JPCPtr->J+1.)*heliAmpCompl;
    _currentParamPreFacMagExpi[*itlamlam] = currentHeliAmp;
    if(!_isWeakDecay){
      std::vector< std::shared_ptr<const JPClamlam> >& currentLPClamlamSymVec=_JPClamlamSymMap.at(*itlamlam);
      std::vector< std::shared_ptr<const JPClamlam> >::iterator itLamLamSym;
      for (itLamLamSym=currentLPClamlamSymVec.begin(); itLamLamSym!=currentLPClamlamSymVec.end(); ++itLamLamSym){
	_currentParamPreFacMagExpi[*itLamLamSym]=(*itLamLamSym)->parityFactor*currentHeliAmp;
      }
    }
  }
   _absDyn->updateFitParams(fitPar);

   if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(fitPar);
   if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(fitPar);
}

void HeliDecAmps::printCurrentAmpParams(Spin& lamX, Spin& lamFs){
  InfoMsg << "Amp name: " << name() << endmsg;
  std::map< std::shared_ptr<const JPClamlam>, complex<double>, pawian::Collection::SharedPtrLess >::iterator it;
  for(it=_currentParamPreFacMagExpi.begin(); it!=_currentParamPreFacMagExpi.end(); ++it){
    Spin lambda1(it->first->lam1);
    if(_enabledlamFsDaughter1 && lamFs!=lambda1) continue;
    Spin lambda2(it->first->lam2);
    if(_enabledlamFsDaughter2 && lamFs!=lambda2) continue;
    Spin lambda(lambda1-lambda2);
    if( fabs(lambda) > it->first->J) continue;

    if (ErrLogger::instance().level()>logging::log_level::INFO) it->first->print(std::cout);
    InfoMsg << " it->first->parityFactor: " << it->first->parityFactor << endmsg;
    InfoMsg << "currentParamPreFacMagExpi: " 
	    << it->second << endmsg;
    }
}


complex<double> HeliDecAmps::XdecPartAmp(const Spin& lamX, Spin& lamDec, short fixDaughterNr, 
					 EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){
  complex<double> result(0.,0.);

  bool lamFs_daughter1=false;
  if( _daughter1IsStable && _Jdaughter1>0) lamFs_daughter1=true;

  bool lamFs_daughter2=false;
  if( _daughter2IsStable && _Jdaughter2>0) lamFs_daughter2=true;

  std::map< std::shared_ptr<const JPClamlam>, complex<double>, pawian::Collection::SharedPtrLess >::iterator it;
  for(it=_currentParamPreFacMagExpi.begin(); it!=_currentParamPreFacMagExpi.end(); ++it){

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

    Id3StringType IdJLamXLam12=FunctionUtils::spin3Index(_J, lamX, lambda);

    complex<double> amp = it->second *                                  
      conj(theData->WignerDIdId3.at(_decay->wigDWigDRefId()).at(IdJLamXLam12));
    result+=amp;
  }

  return result;
}


complex<double> HeliDecAmps::XdecAmp(const Spin& lamX, EvtData* theData, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);
  if( fabs(lamX) > _JPCPtr->J) return result;
  short currentSpinIndex=FunctionUtils::spin1IdIndex(_projId, lamX);

  if (!_recalculate){
    return _cachedAmpIdMap.at(theData->evtNo).at(_absDyn->grandMaId(grandmaAmp)).at(currentSpinIndex);
  }

  std::map< std::shared_ptr<const JPClamlam>, complex<double>, pawian::Collection::SharedPtrLess >::iterator it;
  for(it=_currentParamPreFacMagExpi.begin(); it!=_currentParamPreFacMagExpi.end(); ++it){

    Spin lambda1(it->first->lam1);
    if(_daughter1IsStable && _lam1MinProj!=lambda1) continue;
    Spin lambda2(it->first->lam2);
    if(_daughter2IsStable && _lam2MinProj!=lambda2) continue;
 
    Spin lambda(lambda1-lambda2);
    if( fabs(lambda) > it->first->J) continue;
    
      //      InfoMsg << "_J: " << _J << " lamX: " << lamX << " lambda1: " << lambda1 << " lambda2: " << lambda2
      //	      << " parityFactor: " << _parityFactor << " _currentParamPreFacMagExpi.at(it->first): "
      //	      << _currentParamPreFacMagExpi.at(it->first) << endmsg;
      unsigned int IdJLamXLam12=FunctionUtils::spin3Index(_J, lamX, lambda);

      result+=it->second
        * conj( theData->WignerDIdId3.at(_decay->wigDWigDRefId()).at(IdJLamXLam12) )
	* daughterAmp(lambda1, lambda2, theData);
  }

  if (_absDyn->isLdependent()) result*=_cachedDynL;
  else result*=_cachedDynIdMap.at(_absDyn->grandMaId(grandmaAmp));

      if(result.real()!=result.real()){
        Alert << "result:\t" << result << endmsg;
        exit(0);
      }
      
      if ( _cacheAmps){
	_cachedAmpIdMap[theData->evtNo][_absDyn->grandMaId(grandmaAmp)][currentSpinIndex]=result;
      }

  
  return result;
}

void HeliDecAmps::calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp){
  if(!_recalculate) return; 

  if(!_absDyn->isLdependent()){
    AbsXdecAmp::calcDynamics(theData, grandmaAmp);
    //    return;
  }
  else _cachedDynL=_absDyn->eval(theData, grandmaAmp, absDec()->orbMomMin());

 if(!_daughter1IsStable) _decAmpDaughter1->calcDynamics(theData, this);
 if(!_daughter2IsStable) _decAmpDaughter2->calcDynamics(theData, this);
 return;
}
