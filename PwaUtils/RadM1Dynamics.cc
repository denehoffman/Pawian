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
#include "Utils/IdStringMapRegistry.hh"

RadM1Dynamics::RadM1Dynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::vector<Particle*>& fsParticlesDaughter1, std::vector<Particle*>& fsParticlesDaughter2, double massB0) :
  AbsDynamics(name, fsParticles, mother)
  ,_fsParticlesDaughter1(fsParticlesDaughter1)
  ,_fsParticlesDaughter2(fsParticlesDaughter2)
  ,_dynMassKeyDaughter1(_dynKey+FunctionUtils::particleListName(fsParticlesDaughter1))
  ,_dynMassKeyDaughter2(_dynKey+FunctionUtils::particleListName(fsParticlesDaughter2))
  ,_dynEgammaCMmotherKey(_dynKey+FunctionUtils::particleListName(fsParticlesDaughter1)+"Gamma")
  ,_dynMassIdDaughter1(IdStringMapRegistry::instance()->keyStringId(_keyForMassList, _dynMassKeyDaughter1)) 
  ,_dynMassIdDaughter2(IdStringMapRegistry::instance()->keyStringId(_keyForMassList, _dynMassKeyDaughter2))
  ,_dynEgammaCMmotherId(IdStringMapRegistry::instance()->keyStringId(_keyForMassList, _dynEgammaCMmotherKey))  
  ,_isP1Gamma(true)
  ,_massB0(massB0)
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
    return _cachedMap.at(evtNo);
  }
  
  // Which Daughter particle is the radiative photon?
  double massB = theData->DoubleMassId.at(_dynMassIdDaughter2); // DUMMY VALUE - needs to be mass of non-gamma dacay particle, e.g. eta_c
  if (!_isP1Gamma) massB = theData->DoubleMassId.at(_dynMassIdDaughter1);

  //  double currentMassB = 1.; // DUMMY VALUE - needs to be current mass (fit parameter!) of non-gamma decay particle, e.g. eta_c
  double Egamma = theData->DoubleMassId.at(_dynEgammaCMmotherId); // how to access Egamma?

  complex<double> result(1.,0.);
  result=RadMultipoleFormFactor::PureM1(theData->DoubleMassId.at(_dynId), massB, _massB0, Egamma);
  if ( _cacheAmps) _cachedMap[evtNo]=result;

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
    std::string currentName=(*it)->name();
    mass4VecD1+= theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(currentName));
  }

  double sqrMassDaughter1=mass4VecD1.Mass2();
  if (sqrMassDaughter1>0.) theData->DoubleMassId[_dynMassIdDaughter1]=mass4VecD1.Mass();
  else if( sqrMassDaughter1 > -1.e-6) theData->DoubleMassId[_dynMassIdDaughter1]=0.;
  else{
    Alert << "mass4VecD1.Mass2() is " << mass4VecD1.Mass2() << " and thus < -1e-6 !!!" 
	  << "\nexit !!!" << endmsg;
    exit(0); 
  }

  Vector4<double> mass4VecD2(0.,0.,0.,0.);
  for (it=_fsParticlesDaughter2.begin(); it !=_fsParticlesDaughter2.end(); ++it){
    std::string currentName=(*it)->name();
    mass4VecD2+= theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(currentName));
  }

  double sqrMassDaughter2=mass4VecD2.Mass2();
  if (sqrMassDaughter2>0.) theData->DoubleMassId[_dynMassIdDaughter2]=mass4VecD2.Mass();
  else if( sqrMassDaughter2 > -1.e-6) theData->DoubleMassId[_dynMassIdDaughter2]=0.;
  else{
    Alert << "mass4VecD2.Mass2() is " << mass4VecD2.Mass2() << " and thus < -1e-6 !!!" 
	  << "\nexit !!!" << endmsg;
    exit(0); 
  }

  Vector4<double> mother4Vec(0.,0.,0.,0.);
  mother4Vec=mass4VecD1+mass4VecD2;

  Vector4<double> photonCMmother4Vec(0.,0.,0.,0.);
  if( mass4VecD1.Mass2() < 1.e-12){  //this is the photon work around
    photonCMmother4Vec=mass4VecD1;
  }    
  else if ( mass4VecD2.Mass2() < 1.e-12){  //this is the photon work around 
    _dynEgammaCMmotherKey=_dynMassKeyDaughter2+"Gamma";
    _isP1Gamma=false;
    photonCMmother4Vec=mass4VecD2;
  }
  else{
    Alert << "neither particle 1 with mass " << mass4VecD1.Mass() 
	  << "\tnor particle 2 with mass " << mass4VecD2.Mass()
	  << "\tcan be assigned to the photon!!!" << endmsg;
    exit(1); 
  }
  photonCMmother4Vec.Boost(mother4Vec);
  theData->DoubleMassId[_dynEgammaCMmotherId]=photonCMmother4Vec.E();
}
