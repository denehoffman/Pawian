//************************************************************************//
//									  //
//  Copyright 2019 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// LinearDynamics class definition file. -*- C++ -*-
// Copyright 2019 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/LinearDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "FitParams/AbsPawianParameters.hh"

LinearDynamics::LinearDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, double refMass) :
  AbsDynamics(name, fsParticles, mother),
  _refMass(refMass),
  _currentSlopeMag(0.),
  _currentSlopePhase(0.)
{
  _isLdependent=false;
}

LinearDynamics::~LinearDynamics()
{
}

complex<double> LinearDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  //  int evtNo=theData->evtNo;
  if (!_recalculate){
    return _cachedMap[theData->evtNo];
  }

  double fSlopeMag=-_currentSlopeMag*_refMass+_currentSlopeMag*theData->DoubleMassId.at(_dynId);
  double fSlopePhase=-_currentSlopePhase*_refMass+_currentSlopePhase*theData->DoubleMassId.at(_dynId);
  complex<double> result(complex<double>(1.,0) + std::polar(fSlopeMag, fSlopePhase));
  
  if ( _cacheAmps){
     _cachedMap[theData->evtNo]=result;
  }  

  return result;
}

void  LinearDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  //fill mass
  std::string slopeMagName=_massKey+"SlopeMag";
  fitPar->Add(slopeMagName, 0., 0.001);
  //  fitPar->SetLimits(massName, minMass, maxMass);

  //fill width
  std::string slopePhaseName=_massKey+"SlopePhase";
  fitPar->Add(slopePhaseName, 0., 0.001);
  //  fitPar->SetLimits(widthName, minWidth, maxWidth);
}

void LinearDynamics::fillParamNameList(){
  _paramNameList.clear();

  std::string slopeMagName=_massKey+"SlopeMag";
  _paramNameList.push_back(slopeMagName);

  std::string slopePhaseName=_massKey+"SlopePhase";
  _paramNameList.push_back(slopePhaseName);

}

void LinearDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  std::string slopeMagName=_massKey+"SlopeMag";
  _currentSlopeMag=fitPar->Value(slopeMagName);

  std::string slopePhaseName=_massKey+"SlopePhase";
  _currentSlopePhase=fitPar->Value(slopePhaseName);
}

void LinearDynamics::setMassKey(std::string& theMassKey){
  _massKey=theMassKey;
}
