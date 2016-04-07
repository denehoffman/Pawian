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

// AbsDynamics class definition file. -*- C++ -*-
// Copyright 20123Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/AbsDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/FunctionUtils.hh"
#include "Utils/IdStringMapRegistry.hh"

AbsDynamics::AbsDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsParamHandler()
  ,_name(name)
  ,_massKey("defaultMassKey")
  ,_keyForMassList("Mass")
  ,_fsParticles(fsParticles)
  ,_mother(mother)
  ,_dynKey(FunctionUtils::particleListName(fsParticles))
  ,_dynId(IdStringMapRegistry::instance()->keyStringId(_keyForMassList, _dynKey))
  ,_grandmaKey("default")
  ,_grandmaId(IdStringMapRegistry::instance()->keyStringId("grandMaAndMassParKey", _grandmaKey))
  ,_isLdependent(true)
{
  if(0!=mother) _massKey=mother->name(); 
}

AbsDynamics::~AbsDynamics()
{
}

void AbsDynamics::cacheAmplitudes(){
  _cacheAmps=true;
}

void AbsDynamics::fillMasses(EvtData* theData){

  Vector4<double> mass4Vec(0.,0.,0.,0.);
  std::vector<Particle*>::iterator it;
  for (it=_fsParticles.begin(); it != _fsParticles.end(); ++it){
    std::string currentName=(*it)->name();
    mass4Vec+=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(currentName));
  }

  theData->DoubleMassId[_dynId]=mass4Vec.Mass();
}

