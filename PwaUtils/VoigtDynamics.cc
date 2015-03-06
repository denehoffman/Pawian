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

void  VoigtDynamics::getDefaultParams(fitParCol& fitVal, fitParCol& fitErr){
    fitVal.Masses[_massKey]=_mother->mass();
    fitErr.Masses[_massKey]=0.03;
    fitVal.Widths[_massKey]=_mother->width();
    fitErr.Widths[_massKey]=0.2*_mother->width();
    fitVal.Widths[_massSigmaKey]=0.01;
    fitErr.Widths[_massSigmaKey]=0.4*0.01;
}

bool VoigtDynamics::checkRecalculation(fitParCol& theParamVal){
  _recalculate=false;

  double mass=theParamVal.Masses[_massKey];
  if (!CheckDoubleEquality(mass, _currentMass0)){
    _recalculate=true;
  }
  double width=theParamVal.Widths[_massKey];
  if (!CheckDoubleEquality(width, _currentWidth)){
    _recalculate=true;
  }
  double sigma=theParamVal.Widths[_massSigmaKey];
  if (!CheckDoubleEquality(sigma, _currentSigma)){
    _recalculate=true;
  }
  return _recalculate;
}

void VoigtDynamics::updateFitParams(fitParCol& theParamVal){
  _currentMass0=theParamVal.Masses[_massKey];
  _currentWidth=theParamVal.Widths[_massKey];
  _currentSigma=theParamVal.Widths[_massSigmaKey];
}

void VoigtDynamics::setMassKey(std::string& theMassKey){
  _massKey=theMassKey;
  _massSigmaKey=theMassKey+"Sigma";
}
