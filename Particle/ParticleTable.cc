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

#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"

#include "ErrLogger/ErrLogger.hh"
#include <math.h>

ParticleTable::ParticleTable()
{
}

ParticleTable::~ParticleTable()
{
}

bool ParticleTable::addParticle(Particle* newParticle)
{
  std::string* key = new std::string(newParticle->name());
  particles.insert(pdt_value_type(*key, newParticle));
  return true;
}

Particle* ParticleTable::particle(std::string name)
{
  std::map<const std::string, Particle*>::const_iterator iter = particles.find(name);
  if (iter != particles.end())
    return iter->second;
  else
    return 0;
}

Particle* ParticleTable::particle(double mass)
{
  std::map<const std::string, Particle*>::const_iterator iter = particles.begin();

  while (iter != particles.end() && fabs((iter->second->mass() - mass) / mass) < 0.02)
    ++iter;
  
  if (iter != particles.end())
    return iter->second;
  else
    return 0;
}


bool ParticleTable::clone(std::string newName, std::string oldName)
{
  std::map<const std::string, Particle*>::const_iterator iter = particles.find(newName);
  if (iter != particles.end()) {
    Warning << "ParticleTable: can not clone to already existing particle " << newName << endmsg;
    return false; // cannot clone to an existing particle
  }
  
  Particle* oldParticle = particle(oldName);
  if (0 == oldParticle) {
    Warning << "ParticleTable: cannot clone, " << oldName << " does not exist!" << endmsg;
    return false;
  }

  ParticleData* newData = new ParticleData(oldParticle->data());
  newData->name = newName;
  addParticle(new Particle(*newData));

  return 1; // success
}

bool ParticleTable::modifyMass(std::string name, double newMass)
{
  Particle* thisParticle = particle(name);
  if (0 != thisParticle) {
    thisParticle->data()->mass.mean(newMass);
    return true; // success
  }
  Warning << "ParticleTable::modifyMass: " << name << " not found!" << endmsg;
  return false; // error
}

bool ParticleTable::modifyWidth(std::string name, double newWidth)
{
  Particle* thisParticle = particle(name);
  if (0 != thisParticle) {
    thisParticle->data()->width.mean(newWidth);
    return true; // success
  }
  Warning << "ParticleTable::modifyWidth: " << name << " not found!" << endmsg;
  return false; // error
}


void ParticleTable::print(std::ostream& o) const
{
  std::map<const std::string, Particle*>::const_iterator iter = particles.begin();
  while (iter != particles.end()) {
    iter->second->print(o);
    ++iter;
  }
  return;
}
