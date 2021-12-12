//************************************************************************//
//									  //
//  Copyright 2021 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// VoigtBlattWRelDynamics class definition file. -*- C++ -*-
// Copyright 20213Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/VoigtBlattWRelDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/BreitWignerFunction.hh"
#include "Utils/IdStringMapRegistry.hh"

VoigtBlattWRelDynamics::VoigtBlattWRelDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::vector<Particle*>& fsParticlesDaughter1, std::vector<Particle*>& fsParticlesDaughter2, double qR) :
  BreitWignerBlattWRelDynamics(name, fsParticles, mother, fsParticlesDaughter1, fsParticlesDaughter2, qR)
  ,_massSigmaKey("defaultMassSigmaKey")
  ,_currentSigma(1.)
{
  _isLdependent=true;
}

VoigtBlattWRelDynamics::~VoigtBlattWRelDynamics()
{
}

complex<double> VoigtBlattWRelDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

  int orbMom(OrbMom);

  if (!_recalculate){
    return _cachedLMap.at(theData->evtNo).at(orbMom);
  }
  //the convolution with a Gassian must be put in here
  complex<double> result(std::abs(BreitWignerFunction::BlattWRel(orbMom, theData->DoubleMassId.at(_dynId), _currentMass, _currentWidth, theData->DoubleMassId.at(_dynMassIdDaughter1), theData->DoubleMassId.at(_dynMassIdDaughter2), _qR)),0.);

  if ( _cacheAmps){
     theMutex.lock();
     _cachedLMap[theData->evtNo][orbMom]=result;
     theMutex.unlock();
  }  

  return result;
}

void  VoigtBlattWRelDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  BreitWignerBlattWRelDynamics::fillDefaultParams(fitPar);
    //fill sigma width
  fitPar->Add(_massSigmaKey, 0.01, 0.4*0.01);
  fitPar->SetLimits(_massSigmaKey, 0., 0.06);
}

void VoigtBlattWRelDynamics::fillParamNameList(){
  _paramNameList.clear();
  BreitWignerBlattWRelDynamics::fillParamNameList();

  //fill sigma width
  _paramNameList.push_back(_massSigmaKey);
}

void VoigtBlattWRelDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  BreitWignerBlattWRelDynamics::updateFitParams(fitPar);
  _currentSigma=fitPar->Value(_massSigmaKey);
}

void VoigtBlattWRelDynamics::setMassKey(std::string& theMassKey){
  BreitWignerBlattWRelDynamics::setMassKey(theMassKey);
  _massSigmaKey=theMassKey+"Sigma";
}
