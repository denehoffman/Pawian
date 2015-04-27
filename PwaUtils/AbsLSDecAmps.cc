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

// AbsLSDecAmps class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/AbsLSDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"
#include "FitParams/AbsPawianParameters.hh"

AbsLSDecAmps::AbsLSDecAmps(std::shared_ptr<IsobarLSDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_LSs(theDec->LSAmps())
  ,_factorMag(1.)
  ,_lam1Min(-_Jdaughter1)
  ,_lam1Max(_Jdaughter1)
  ,_lam2Min(-_Jdaughter2)
  ,_lam2Max(_Jdaughter2)
{
  if(_LSs.size()>0) _factorMag=1./sqrt(_LSs.size());
  Particle* daughter1=_decay->daughter1Part();
  Particle* daughter2=_decay->daughter2Part();
  _parityFactor=daughter1->theParity()*daughter2->theParity()*pow(-1,_JPCPtr->J-daughter1->J()-daughter2->J());
  Info << "_parityFactor=\t" << _parityFactor << endmsg;
  fillCgPreFactor();
}

AbsLSDecAmps::AbsLSDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_lam1Min(-_Jdaughter1)
  ,_lam1Max(_Jdaughter1)
  ,_lam2Min(-_Jdaughter2)
  ,_lam2Max(_Jdaughter2)
{
  Particle* daughter1=_decay->daughter1Part();
  Particle* daughter2=_decay->daughter2Part();
  _parityFactor=daughter1->theParity()*daughter2->theParity()*pow(-1,_JPCPtr->J-daughter1->J()-daughter2->J());
  Info << "_parityFactor=\t" << _parityFactor << endmsg;
  fillCgPreFactor();
}

AbsLSDecAmps::~AbsLSDecAmps()
{
}

void  AbsLSDecAmps::getDefaultParams(fitParCol& fitVal, fitParCol& fitErr){

  std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
  std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
  std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
  std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    currentMagValMap[*itLS]=_factorMag;
    currentPhiValMap[*itLS]=0.;
    currentMagErrMap[*itLS]=_factorMag;
    currentPhiErrMap[*itLS]=0.3;
  }

  fitVal.MagsLS[_key]=currentMagValMap;
  fitVal.PhisLS[_key]=currentPhiValMap;
  fitErr.MagsLS[_key]=currentMagErrMap;
  fitErr.PhisLS[_key]=currentPhiErrMap;

  _absDyn->getDefaultParams(fitVal, fitErr);


  if(!_daughter1IsStable) _decAmpDaughter1->getDefaultParams(fitVal, fitErr);
  if(!_daughter2IsStable) _decAmpDaughter2->getDefaultParams(fitVal, fitErr);
}

void  AbsLSDecAmps::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){

  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    //fill magnitude
    std::string magName=(*itLS)->name()+_key+"Mag";
    double valMag=_factorMag;
    double errMag=_factorMag/2.;
    double minMag=0.;
    double maxMag=_factorMag+30.*errMag;

    fitPar->Add(magName, valMag, errMag);
    fitPar->SetLimits(magName, minMag, maxMag);

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


void AbsLSDecAmps::print(std::ostream& os) const{
  return; //dummy
}


bool AbsLSDecAmps::checkRecalculation(fitParCol& theParamVal){
  _recalculate=false;

   if(_absDyn->checkRecalculation(theParamVal)) _recalculate=true;

   if(!_daughter1IsStable) {
     if(_decAmpDaughter1->checkRecalculation(theParamVal)) _recalculate=true;
   }
   if(!_daughter2IsStable){
     if(_decAmpDaughter2->checkRecalculation(theParamVal)) _recalculate=true;
   }

   if(!_recalculate){
     std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& magMap=theParamVal.MagsLS[_key];
     std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& phiMap=theParamVal.PhisLS[_key];

     std::vector< std::shared_ptr<const LScomb> >::iterator it;
     for (it=_LSs.begin(); it!=_LSs.end(); ++it){
       double theMag=magMap[*it];
       double thePhi=phiMap[*it];

       if(!CheckDoubleEquality(theMag, _currentParamMags[*it])){
	 _recalculate=true;
	 return _recalculate;
       }
       if(!CheckDoubleEquality(thePhi, _currentParamPhis[*it])){
	 _recalculate=true;
	 return _recalculate;
       }
     }
   }

   return _recalculate;
}


void  AbsLSDecAmps::updateFitParams(fitParCol& theParamVal){
   std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& magMap=theParamVal.MagsLS[_key];
   std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& phiMap=theParamVal.PhisLS[_key];

   std::vector< std::shared_ptr<const LScomb> >::iterator it;
   for (it=_LSs.begin(); it!=_LSs.end(); ++it){
     double theMag=magMap[*it];
     double thePhi=phiMap[*it];
     _currentParamMags[*it]=theMag;
     _currentParamPhis[*it]=thePhi;
     complex<double> expi(cos(thePhi), sin(thePhi));
     _currentParamMagExpi[*it]=theMag*expi;
     _currentParamPreFacMagExpi[*it]=_preFactor*_isospinCG*theMag*expi;
   }

   _absDyn->updateFitParams(theParamVal);

  if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(theParamVal);
  if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(theParamVal);

}

void  AbsLSDecAmps::fillCgPreFactor(){

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
