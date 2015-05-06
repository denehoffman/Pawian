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
#include "FitParams/AbsPawianParameters.hh"

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

  return result;
}




complex<double> TensorOmegaTo3PiDecAmps::XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);

  int evtNo=theData->evtNo;
  Id2StringType currentSpinIndex=FunctionUtils::spin2Index(lamX,lamFs);
  
  if ( _cacheAmps && !_recalculate){
    result=_cachedAmpMap.at(evtNo).at(_absDyn->grandMaKey(grandmaAmp)).at(currentSpinIndex);
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
     _cachedAmpMap[evtNo][_absDyn->grandMaKey(grandmaAmp)][currentSpinIndex]=result;
     theMutex.unlock();
  }

  result*=_absDyn->eval(theData, grandmaAmp);
  //  Info <<"TensorOmegaTo3PiDecAmps result: " << result << endmsg; 
  return result;
}

void TensorOmegaTo3PiDecAmps::print(std::ostream& os) const{
  return; //dummy
}


void TensorOmegaTo3PiDecAmps::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){

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

void TensorOmegaTo3PiDecAmps::fillParamNameList(){
  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    std::string magName=(*itLS)->name()+_key+"Mag";
    _paramNameList.push_back(magName);

    std::string phiName=(*itLS)->name()+_key+"Phi";
    _paramNameList.push_back(phiName);
  }
}

void TensorOmegaTo3PiDecAmps::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){

  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    //fill magnitude
    std::string magName=(*itLS)->name()+_key+"Mag";
    double theMag= fabs(fitPar->Value(magName));
    _currentParamMags[*itLS]=theMag;

    std::string phiName=(*itLS)->name()+_key+"Phi";
    double thePhi= fitPar->Value(phiName);
    _currentParamPhis[*itLS]=thePhi;
  }

  _absDyn->updateFitParams(fitPar);

  if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(fitPar);
  if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(fitPar);
}

