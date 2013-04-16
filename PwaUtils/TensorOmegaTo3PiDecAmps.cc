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

// TensorOmegaTo3PiDecAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/TensorOmegaTo3PiDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/OmegaTo3PiTensorDecay.hh"
#include "Particle/Particle.hh"


TensorOmegaTo3PiDecAmps::TensorOmegaTo3PiDecAmps(boost::shared_ptr<OmegaTo3PiTensorDecay> theDec) :
  AbsXdecAmp(theDec)
  ,_JPCLSs(theDec->JPCLSAmps())
  ,_factorMag(1.)
{
  if(_JPCLSs.size()>0) _factorMag=1./sqrt(_JPCLSs.size());
  _daughter1=_decay->daughter1Part();
  _daughter2=_decay->daughter2Part();
  _daughter3=theDec->daughter3Part();
}

TensorOmegaTo3PiDecAmps::~TensorOmegaTo3PiDecAmps()
{
}


complex<double> TensorOmegaTo3PiDecAmps::XdecPartAmp(Spin lamX, Spin lamDec, short fixDaughterNr, EvtData* theData, Spin lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);
  // std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
  // for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
  //   if( fabs(lamX) > (*it)->J ) continue;
  //   double theMag=_currentParamMags[*it];
  //   double thePhi=_currentParamPhis[*it];
  //   complex<double> expi(cos(thePhi), sin(thePhi));

  //       complex<double> amp = theMag*expi*sqrt(2*(*it)->L+1)
  // 	  *conj( theData->WignerDsString[_wignerDKey][(*it)->J][lamX][0]);
  //       result+=amp;
  // }

  return result;
}




complex<double> TensorOmegaTo3PiDecAmps::XdecAmp(Spin lamX, EvtData* theData, Spin lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.); 
  
  int evtNo=theData->evtNo;

  if ( _cacheAmps && !_recalculate){
    result=_cachedAmpMap[evtNo][lamX][lamFs];
    return result;
  }

 
  // Vector4<double> part1_4Vec=theData->FourVecsString[_daughter1->name()];
  // Vector4<double> part2_4Vec=theData->FourVecsString[_daughter2->name()];
  // Vector4<double> part3_4Vec=theData->FourVecsString[_daughter3->name()];
 
  // Vector4<double> P_2part_4Vec=part1_4Vec+part2_4Vec+part3_4Vec;
  // LeviCivitaTensor eps; 
  // PolVector omega; // spin-1 particle is the default constructor
  // omega.SetP4(P_2part_4Vec, P_2part_4Vec.M());
  
  std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
  for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
    if( fabs(lamX) > (*it)->J ) continue;
    double theMag=_currentParamMags[*it];
    double thePhi=_currentParamPhis[*it];
    complex<double> expi(cos(thePhi), sin(thePhi));
    //    Tensor<complex<double> >  ampTensor(0);
    //    ampTensor = eps|(part1_4Vec%part2_4Vec%part3_4Vec%omega(lamX));   
    complex<double> amp = theMag*expi*theData->ComplexDoubleString["omegTensor"][_JPCPtr->J][lamX];

    //    Info << "amp:\t" << amp << endmsg;

    result+=amp;
  }


  if ( _cacheAmps){
     theMutex.lock();
     _cachedAmpMap[evtNo][lamX][lamFs]=result;
     theMutex.unlock();
  }

  result*=_absDyn->eval(theData, grandmaAmp);
  return result;
}

complex<double> TensorOmegaTo3PiDecAmps::lsLoop(Spin lamX, EvtData* theData, Spin lam1Min, Spin lam1Max, Spin lam2Min, Spin lam2Max, bool withDecs, Spin lamFs ){
  complex<double> result(0.,0.);
  // std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
  // for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){

  //   double theMag=_currentParamMags[*it];
  //   double thePhi=_currentParamPhis[*it];
  //   complex<double> expi(cos(thePhi), sin(thePhi));

  //   for(Spin lambda1=lam1Min; lambda1<=lam1Max; ++lambda1){
  //     for(Spin lambda2=lam2Min; lambda2<=lam2Max; ++lambda2){
  //       Spin lambda = lambda1-lambda2;
  //       if( fabs(lambda)>(*it)->J || fabs(lambda)>(*it)->S) continue;
  //       complex<double> amp = theMag*expi*_cgPreFactor[*it][lambda1][lambda2]*conj( theData->WignerDsString[_wignerDKey][(*it)->J][lamX][lambda]);

  //       if(withDecs) amp *=daughterAmp(lambda1, lambda2, theData, lamFs, this);
  //       result+=amp;
  //     }
  //   }
  // }
  return result;
}

void TensorOmegaTo3PiDecAmps::print(std::ostream& os) const{
  return; //dummy
}


void TensorOmegaTo3PiDecAmps::getDefaultParams(fitParams& fitVal, fitParams& fitErr){

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
  for(itLS=_JPCLSs.begin(); itLS!=_JPCLSs.end(); ++itLS){
    currentMagValMap[*itLS]=_factorMag;
    currentPhiValMap[*itLS]=0.;
    currentMagErrMap[*itLS]=_factorMag;
    currentPhiErrMap[*itLS]=0.3;
  }

  fitVal.Mags[_key]=currentMagValMap;
  fitVal.Phis[_key]=currentPhiValMap;
  fitErr.Mags[_key]=currentMagErrMap;
  fitErr.Phis[_key]=currentPhiErrMap;

  _absDyn->getDefaultParams(fitVal, fitErr);


  if(!_daughter1IsStable) _decAmpDaughter1->getDefaultParams(fitVal, fitErr);
  if(!_daughter2IsStable) _decAmpDaughter2->getDefaultParams(fitVal, fitErr);  
}


bool TensorOmegaTo3PiDecAmps::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;

   if(_absDyn->checkRecalculation(theParamVal)) _recalculate=true; 

   if(!_daughter1IsStable) {
     if(_decAmpDaughter1->checkRecalculation(theParamVal)) _recalculate=true;
   }
   if(!_daughter2IsStable){
     if(_decAmpDaughter2->checkRecalculation(theParamVal)) _recalculate=true;
   }

   if(!_recalculate){
     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& magMap=theParamVal.Mags[_key];
     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& phiMap=theParamVal.Phis[_key];
     
     std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
     for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
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


void TensorOmegaTo3PiDecAmps::updateFitParams(fitParams& theParamVal){
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& magMap=theParamVal.Mags[_key];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& phiMap=theParamVal.Phis[_key];

   std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
   for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
     double theMag=magMap[*it];
     double thePhi=phiMap[*it];
     _currentParamMags[*it]=theMag;
     _currentParamPhis[*it]=thePhi;
   }

   _absDyn->updateFitParams(theParamVal);

  if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(theParamVal);
  if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(theParamVal);

}




