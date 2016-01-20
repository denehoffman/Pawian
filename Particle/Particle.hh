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

#include <string>
#include <iostream>
#include "Particle/Parity.hh"
#include "Utils/ErrValue.hh"

enum ParticleType {parton, lepton, boson, meson, baryon, nucleus, special};
enum DynFunctionType {BW, relBW, relBWBlattWK, undefined};

struct ParticleData {
  ParticleType type;
  std::string name;
  std::string texName;
  int charge;
  ErrValue mass;
  ErrValue width;
  int twoJ;
  int twoIso;
  int twoIsoThree;
  int strange;
  int charm;
  int theParity;
  Parity parity;
  int theCParity;
  Parity chargeParity;
  Parity gParity;
  int theGParity;
  DynFunctionType dynamicFunction;
  ParticleData(){}
  ParticleData(ParticleData* other):
    type(other->type),
    name(other->name),
    texName(other->texName),
    charge(other->charge),
    mass(other->mass),
    width(other->width),
    twoJ(other->twoJ),
    twoIso(other->twoIso),
    twoIsoThree(other->twoIsoThree),
    strange(other->strange),
    charm(other->charm), 
    theParity(other->theParity),
    parity(other->parity),
    theCParity(other->theCParity),
    gParity(other->gParity),
    theGParity(other->theGParity),
    dynamicFunction(other->dynamicFunction)
  { }
};


class Particle
{
public:

  Particle();
  Particle(const ParticleData& data);
  Particle(const Particle& other);
  ~Particle();

  const std::string& name();
  const std::string& texName();
  ParticleType type();
  int charge();
  const ErrValue& massErr() const;
  const ErrValue& widthErr() const;
  double mass();
  double width();
  int twoJ();
  double J();
  int theParity();
  Parity& parity();
  int theCParity();
  Parity& chargeParity();
  Parity& gParity();
  int theGParity();
  int twoIso();
  int twoIso3();
  int strange();
  int charm();
  DynFunctionType dynFctType();
  ParticleData* data() const;
  bool sameName(Particle& compare);
  bool operator<(Particle& compare);
  bool operator==(Particle& compare);
  void print(std::ostream& out);

private:
  ParticleData* pdata;

friend std::ostream &operator<<(std::ostream &o, Particle &p);
};
