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
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/OmegaTo3PiTensorDecay.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

TensorOmegaTo3PiDecAmps::TensorOmegaTo3PiDecAmps(std::shared_ptr<OmegaTo3PiTensorDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_LSs(theDec->LSAmps())
  ,_factorMag(1.)
{
  if(_LSs.size()>0) _factorMag=1./sqrt(_LSs.size());
  _daughter1=_decay->daughter1Part();
  _daughter2=_decay->daughter2Part();
  _daughter3=theDec->daughter3Part();
}

TensorOmegaTo3PiDecAmps::~TensorOmegaTo3PiDecAmps()
{
}


complex<double> TensorOmegaTo3PiDecAmps::XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);
  // std::vector< std::shared_ptr<const JPCLS> >::iterator it;
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




complex<double> TensorOmegaTo3PiDecAmps::XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);

  int evtNo=theData->evtNo;
  unsigned short currentSpinIndex=lamX.ToIndex()*100+lamFs.ToIndex();

  if ( _cacheAmps && !_recalculate){
    result=_cachedAmpShortMap.at(evtNo).at(currentSpinIndex);
    //    result=_cachedAmpMap[evtNo][lamX][lamFs];
    return result;
  }

  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=_LSs.begin(); it!=_LSs.end(); ++it){
    //    if( fabs(lamX) > _JPCPtr->J ) continue;
    double theMag=_currentParamMags[*it];
    double thePhi=_currentParamPhis[*it];
    complex<double> expi(cos(thePhi), sin(thePhi));

    complex<double> amp = theMag*expi*theData->ComplexDoubleString[_wignerDKey][_JPCPtr->J][lamX];

    result+=amp;
  }

  result*=100.;
  if ( _cacheAmps){
     theMutex.lock();
     _cachedAmpShortMap[evtNo][currentSpinIndex]=result;
     //     _cachedAmpMap[evtNo][lamX][lamFs]=result;
     theMutex.unlock();
  }

  result*=_absDyn->eval(theData, grandmaAmp);
  //  Info <<"TensorOmegaTo3PiDecAmps result: " << result << endmsg; 
  return result;
}

void TensorOmegaTo3PiDecAmps::print(std::ostream& os) const{
  return; //dummy
}


void TensorOmegaTo3PiDecAmps::getDefaultParams(fitParams& fitVal, fitParams& fitErr){

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


void TensorOmegaTo3PiDecAmps::updateFitParams(fitParams& theParamVal){
   std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& magMap=theParamVal.MagsLS[_key];
   std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& phiMap=theParamVal.PhisLS[_key];

   std::vector< std::shared_ptr<const LScomb> >::iterator it;
   for (it=_LSs.begin(); it!=_LSs.end(); ++it){
     double theMag=magMap[*it];
     double thePhi=phiMap[*it];
     _currentParamMags[*it]=theMag;
     _currentParamPhis[*it]=thePhi;
   }

   _absDyn->updateFitParams(theParamVal);

  if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(theParamVal);
  if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(theParamVal);

}




