//************************************************************************//
//									  //
//  Copyright 2025 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// SigmaParamDynamics class definition file. -*- C++ -*-
// Copyright 2025 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/SigmaParamDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/SigmaParameterization.hh"

SigmaParamDynamics::SigmaParamDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsDynamics(name, fsParticles, mother)
  ,_currentStrength(1.)
  ,_sigmaParamDyn(std::shared_ptr<SigmaParameterization>(new SigmaParameterization))
{
  _isLdependent=false;
}

SigmaParamDynamics::~SigmaParamDynamics()
{
}

complex<double> SigmaParamDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  if ( !_recalculate){
    return _cachedMap[theData->evtNo];
  }

  complex<double> result=_sigmaParamDyn->calc(theData->DoubleMassId.at(_dynId), _currentStrength);
  if ( _cacheAmps){
     _cachedMap[theData->evtNo]=result;
  }  

  return result;
}


void SigmaParamDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  //fill mass
  std::string strengthName=_massKey+"Strength";
  fitPar->Add(strengthName, _currentStrength, 0.1);
  //  fitPar->SetLimits(massName, minMass, maxMass);
}


void SigmaParamDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  std::string strengthName=_massKey+"Strength";
  _currentStrength=fitPar->Value(strengthName);
}

void SigmaParamDynamics::fillParamNameList(){
  _paramNameList.clear();
  //fill
  std::string strengthName=_massKey+"Strength";
  _paramNameList.push_back(strengthName);
}
