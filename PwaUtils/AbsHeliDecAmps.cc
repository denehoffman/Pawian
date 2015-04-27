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

// AbsHeliDecAmps class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/AbsHeliDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"
#include "FitParams/AbsPawianParameters.hh"

AbsHeliDecAmps::AbsHeliDecAmps(std::shared_ptr<IsobarHeliDecay> theDec, ChannelID channelID) :
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

AbsHeliDecAmps::AbsHeliDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
{
  Particle* daughter1=_decay->daughter1Part();
  Particle* daughter2=_decay->daughter2Part();
  _parityFactor=_JPCPtr->P*daughter1->theParity()*daughter2->theParity()*pow(-1,_JPCPtr->J-daughter1->J()-daughter2->J());
  Info << "_parityFactor=\t" << _parityFactor << endmsg;
}

AbsHeliDecAmps::~AbsHeliDecAmps()
{
}

void  AbsHeliDecAmps::getDefaultParams(fitParCol& fitVal, fitParCol& fitErr){

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

void  AbsHeliDecAmps::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){

  std::vector< std::shared_ptr<const JPClamlam> >::const_iterator itlamlam;
  for(itlamlam=_JPClamlams.begin(); itlamlam!=_JPClamlams.end(); ++itlamlam){
    //fill magnitude
    std::string magName=(*itlamlam)->name()+_key+"Mag";
    double valMag=_factorMag;
    double errMag=_factorMag/2.;
    double minMag=0.;
    double maxMag=_factorMag+30.*errMag;

    fitPar->Add(magName, valMag, errMag);
    fitPar->SetLimits(magName, minMag, maxMag);

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

void AbsHeliDecAmps::print(std::ostream& os) const{
  return; //dummy
}


bool AbsHeliDecAmps::checkRecalculation(fitParCol& theParamVal){
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


void  AbsHeliDecAmps::updateFitParams(fitParCol& theParamVal){
   std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >& magMap=theParamVal.MagLamLams[_key];
   std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >& phiMap=theParamVal.PhiLamLams[_key];

   std::vector< std::shared_ptr<const JPClamlam> >::iterator it;
   for (it=_JPClamlams.begin(); it!=_JPClamlams.end(); ++it){
     double theMag=magMap[*it];
     double thePhi=phiMap[*it];
     _currentParamMagLamLams[*it]=theMag;
     _currentParamPhiLamLams[*it]=thePhi;
     complex<double> expi(cos(thePhi), sin(thePhi));

     _currentParamMagExpi[*it]=theMag*expi;
     _currentParamPreFacMagExpi[*it]=_preFactor*_isospinCG*sqrt(2.*_JPCPtr->J+1.)*theMag*expi;

     std::vector< std::shared_ptr<const JPClamlam> >& currentLPClamlamVec=_JPClamlamSymMap[*it];
     std::vector< std::shared_ptr<const JPClamlam> >::iterator itLamLam;
     for (itLamLam=currentLPClamlamVec.begin(); itLamLam!=currentLPClamlamVec.end(); ++itLamLam){
       _currentParamMagLamLams[*itLamLam]=theMag;
       _currentParamPhiLamLams[*itLamLam]=thePhi;
       _currentParamMagExpi[*itLamLam]=theMag*expi;
       _currentParamPreFacMagExpi[*itLamLam]=_preFactor*_isospinCG*sqrt(2.*_JPCPtr->J+1.)*theMag*expi;
     }
   }


   _absDyn->updateFitParams(theParamVal);

   if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(theParamVal);
   if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(theParamVal);

}
