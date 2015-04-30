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

// BreitWignerDynamics class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/RadM1Dynamics.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/RadMultipoleFormFactor.hh"
#include "ConfigParser/ParserBase.hh"

RadM1Dynamics::RadM1Dynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::vector<Particle*>& fsParticlesDaughter1, std::vector<Particle*>& fsParticlesDaughter2) :
  AbsDynamics(name, fsParticles, mother)
  ,_fsParticlesDaughter1(fsParticlesDaughter1)
  ,_fsParticlesDaughter2(fsParticlesDaughter2)
  ,_dynMassKeyDaughter1(_dynKey+FunctionUtils::particleListName(fsParticlesDaughter1))
  ,_dynMassKeyDaughter2(_dynKey+FunctionUtils::particleListName(fsParticlesDaughter2))
{
  Info << "RadM1Dynamics for " << _name <<endmsg;
  _isLdependent=false;
}

RadM1Dynamics::~RadM1Dynamics()
{
}

complex<double> RadM1Dynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate){
    return _cachedMap[evtNo];
  }
  
  // Which Daughter particle is the radiative photon?
  double massB = 1.; // DUMMY VALUE - needs to be mass of non-gamma dacay particle, e.g. eta_c
  double currentMassB = 1.; // DUMMY VALUE - needs to be current mass (fit parameter!) of non-gamma decay particle, e.g. eta_c
  double Egamma = 1.; // how to access Egamma?

  complex<double> result(1.,0.);
  result=RadMultipoleFormFactor::PureM1(theData->DoubleString.at(_dynKey), massB, currentMassB, Egamma);

  return result;
}

void  RadM1Dynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  return;
}

void RadM1Dynamics::fillParamNameList(){
  return;
}

void RadM1Dynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  return;
}

void RadM1Dynamics::setMassKey(std::string& theMassKey){
//  _massKey=theMassKey;
}

void RadM1Dynamics::fillMasses(EvtData* theData){
  AbsDynamics::fillMasses(theData);

  Vector4<double> mass4VecD1(0.,0.,0.,0.);
  std::vector<Particle*>::iterator it;
  for (it=_fsParticlesDaughter1.begin(); it !=_fsParticlesDaughter1.end(); ++it){
    mass4VecD1+=theData->FourVecsString[(*it)->name()];
  }
  theData->DoubleString[_dynMassKeyDaughter1]=mass4VecD1.Mass();

  Vector4<double> mass4VecD2(0.,0.,0.,0.);
  for (it=_fsParticlesDaughter2.begin(); it !=_fsParticlesDaughter2.end(); ++it){
    mass4VecD2+=theData->FourVecsString[(*it)->name()];
  }
  theData->DoubleString[_dynMassKeyDaughter2]=mass4VecD2.Mass();
}
