//************************************************************************//
//									  //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#pragma once

#include "ErrLogger/ErrLogger.hh"
#include <string>
#include <vector>

class Particle;

class PhpGenDynamics
{
public:
  PhpGenDynamics(std::string& description, std::vector<Particle*>& finalStateParticles);
  ~PhpGenDynamics();

  double eval(double currentMass);
  const double maxWeight() const {return _maxWeight;}
  bool isDecayParticle(Particle* compParticle);
  std::vector<unsigned int> particleIds() {return _iDparticles;}

private:
  std::string _dynName;
  double _mass0;
  double _width0;
  double _massDaughter1;
  double _massDaughter2;
  double _maxWeight;
  std::vector<unsigned int> _iDparticles;
  std::vector<Particle*> _decParticles;
};
