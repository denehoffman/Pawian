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

// VoigtDynamics class definition file. -*- C++ -*-
// Copyright 20123Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/VoigtDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/Voigtian.hh"
#include "FitParams/AbsPawianParameters.hh"

VoigtDynamics::VoigtDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsDynamics(name, fsParticles, mother),
  _massSigmaKey("defaultMassSigmaKey"),
  _voigtPtr(new Voigtian())
{
  if(0!=mother) _massSigmaKey=_massKey+"Sigma";
  _isLdependent=false;
}

VoigtDynamics::~VoigtDynamics()
{
}

complex<double> VoigtDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate){
    return _cachedMap[evtNo];
  }

  complex<double> result(sqrt(_voigtPtr->calc(theData->DoubleString.at(_dynKey), _currentMass0, _currentWidth, _currentSigma)), 0.);
  
  if ( _cacheAmps){
     theMutex.lock();
     _cachedMap[evtNo]=result;
     theMutex.unlock();
  }  

  return result;
}

void  VoigtDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  //fill mass
  std::string massName=_massKey+"Mass";
  double valMass=_mother->mass();
  double errMass=0.03;
  double minMass=valMass-5.*errMass;
  if(minMass<0.) minMass=0.;
  double maxMass=valMass+5.*errMass;
  
  fitPar->Add(massName, valMass, errMass);
  fitPar->SetLimits(massName, minMass, maxMass);

  //fill width
  std::string widthName=_massKey+"Width";
  double valWidth=_mother->width();
  double errWidth=0.2*_mother->width();
  double minWidth=valWidth-5.*errWidth;
  if(minWidth<0.) minWidth=0.;
  double maxWidth=valWidth+5.*errWidth;

  fitPar->Add(widthName, valWidth, errWidth);
  fitPar->SetLimits(widthName, minWidth, maxWidth);

  //fill sigma width
  fitPar->Add(_massSigmaKey, 0.01, 0.4*0.01);
  fitPar->SetLimits(_massSigmaKey, 0., 0.06);

}

void VoigtDynamics::fillParamNameList(){
  _paramNameList.clear();

  std::string massName=_massKey+"Mass";
  _paramNameList.push_back(massName);

  //fill width
  std::string widthName=_massKey+"Width";
  _paramNameList.push_back(widthName);

  //fill sigma width
  _paramNameList.push_back(_massSigmaKey);
}

void VoigtDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  std::string massName=_massKey+"Mass";
  _currentMass0=fitPar->Value(massName);

  std::string widthName=_massKey+"Width";
  _currentWidth=fitPar->Value(widthName);

  _currentSigma=fitPar->Value(_massSigmaKey);
}

void VoigtDynamics::setMassKey(std::string& theMassKey){
  _massKey=theMassKey;
  _massSigmaKey=theMassKey+"Sigma";
}
