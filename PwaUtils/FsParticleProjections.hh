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

// WoDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class Particle;

class FsParticleProjections{

public:
  FsParticleProjections(std::vector<Particle*>& fsParticles);
  virtual ~FsParticleProjections();

  unsigned int fsParticleId(std::string& name);
  std::vector<std::string> finalStateNames() const {return _fsParticleNames;}
  std::vector< std::vector<Spin> > spinProjections() {return _fsSpinProjections;}

protected:

private:
  std::vector<Particle*> _fsParticles;
  std::vector<std::string> _fsParticleNames;
  std::vector<Spin> _fsSpins;
  std::vector< std::vector<Spin> > _fsSpinProjections;
  bool increaseCurrentSpinIndex(std::vector<unsigned int >& currentSpinProjections, std::vector<unsigned int >& spinMaxArrays, unsigned int& currentPIndex);
};
