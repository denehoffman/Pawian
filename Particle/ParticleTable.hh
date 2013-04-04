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

#pragma once

#include <map>
#include <string>

class Particle;

typedef std::pair<const std::string, Particle*> pdt_value_type;

class ParticleTable
{
public:

  ParticleTable();
  ~ParticleTable();

  Particle* particle(std::string name);
  Particle* particle(double mass);

  bool addParticle(Particle* newParticle);
  bool clone(std::string newName, std::string oldName);

  bool modifyMass(std::string name, double newMass);
  bool modifyWidth(std::string name, double newWidth);

  void print(std::ostream& o) const;

private:
  std::map<const std::string, Particle*> particles;

};

