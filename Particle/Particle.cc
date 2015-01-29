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

#include "Particle/Particle.hh"

#include "ErrLogger/ErrLogger.hh"

#include <cstdlib>
#include <iostream>
#include <string>
#include <cstdlib>
#include <math.h>

Particle::Particle(const ParticleData& data)
{
  pdata = const_cast<ParticleData*>(&data);
}

Particle::Particle(const Particle& other)
{
  pdata = other.data();
}

Particle::~Particle()
{
}

const std::string& Particle::name()
{
  if (0 != pdata)
    return pdata->name;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

const std::string& Particle::texName()
{
  if (0 != pdata)
    return pdata->texName;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

ParticleType Particle::type()
{
  if (0 != pdata)
    return pdata->type;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

int Particle::charge()
{
  if (0 != pdata)
    return pdata->charge;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
  return 0;
}

const ErrValue& Particle::massErr() const
{
  if (0 != pdata)
    return pdata->mass;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

const ErrValue& Particle::widthErr() const
{
  if (0 != pdata)
    return pdata->width;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

double Particle::mass()
{
  if (0 != pdata)
    return (pdata->mass).mean();
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
  return 0.;
}

double Particle::width()
{
  if (0 != pdata)
    return (pdata->width).mean();
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
  return 0.;
}

int Particle::twoJ()
{
  if (0 != pdata)
    return pdata->twoJ;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
  return 0;

}

double Particle::J()
{
  if (0 != pdata)
    return double(pdata->twoJ)/2.;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
  return 0.;
}

int Particle::theParity()
{
  if (0 != pdata)
    return pdata->theParity;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

Parity& Particle::parity()
{
  if (0 != pdata)
    return pdata->parity;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

int Particle::theCParity()
{
  if (0 != pdata)
    return pdata->theCParity;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

Parity& Particle::chargeParity()
{
  if (0 != pdata)
    return pdata->chargeParity;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

Parity& Particle::gParity()
{
  if (0 != pdata)
    return pdata->gParity;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

int Particle::theGParity()
{
  if (0 != pdata)
    return pdata->theGParity;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

int Particle::twoIso()
{
  if (0 != pdata)
    return pdata->twoIso;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

int Particle::twoIso3()
{
  if (0 != pdata)
    return pdata->twoIsoThree;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

int Particle::strange()
{
  if (0 != pdata)
    return pdata->strange;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

int Particle::charm()
{
  if (0 != pdata)
    return pdata->charm;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

DynFunctionType Particle::dynFctType()
{
  if (0 != pdata)
    return pdata->dynamicFunction;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}


ParticleData* Particle::data() const
{
  return pdata;
}

void Particle::print(std::ostream& out)
{
  out << name() << "\tmass=" << massErr() << "\twidth=" << widthErr() << "\t3*q=" << charge()
      << "\t2*J=" << twoJ() << "\tP=" << theParity() << "\tp=" << parity().parity() 
      << "\tC=" << theCParity()  << "\tc=" << chargeParity().parity() 
      << "\tG=" << gParity().parity() << "\ttheGParity=" << theGParity() << "\t2*I=" << twoIso() << "\t2*I3=" << twoIso3() 
      << "\tcharm=" << charm() << "\tstrange=" << strange() << std::endl;
  
  return;
}

Particle::Particle()
{
}

bool Particle::sameName(Particle& compare){
  return ( name()==compare.name() );
}

bool Particle::operator<(Particle& compare){
 bool result=false;
 if( type() < compare.type()) result=true;
 else if( charm() < compare.charm()) result=true;
 else if( strange() < compare.strange()) result=true;
 else if( twoIso() < compare.twoIso()) result=true;
 else if( twoIso3() < compare.twoIso3()) result=true;
 else if( J() < compare.J()) result=true;
 else if( theParity() < compare.theParity()) result=true;
 else if( theCParity() < compare.theCParity()) result=true; 
 else if( theGParity() < compare.theGParity()) result=true;
 else if( charge() < compare.charge()) result=true;
 else if ( mass() < compare.mass()) result=true;
 else if ( width() < compare.width()) result=true;
 return result;
}

bool Particle::operator==(Particle& compare){
 bool result=true;
 if( type() != compare.type()) result=false;
 else if( charm() != compare.charm()) result=false;
 else if( strange() != compare.strange()) result=false;
 else if( twoIso() != compare.twoIso()) result=false;
 else if( twoIso3() != compare.twoIso3()) result=false;
 else if( fabs(J()-compare.J()) > 1.e-8) result=false;
 else if( theParity() != compare.theParity()) result=false;
 else if( theCParity() != compare.theCParity()) result=false;
 else if( theGParity() != compare.theGParity()) result=false;
 else if( charge() != compare.charge()) result=false;
 else if ( fabs(mass()-compare.mass()) > 1.e-8) result=false;
 else if ( fabs(width()-compare.width()) > 1.e-8) result=false;
 return result;
}


std::ostream &operator<<(std::ostream &o, Particle &p)
{
  o << p.name() << "\tmass=" << p.massErr() << "\twidth=" << p.widthErr() << "\t3*q=" << p.charge()
    << "\t2*J=" << p.twoJ() << "\tP=" << p.theParity() << "\tp=" << p.parity().parity() 
    << "\tC=" << p.theCParity()  << "\tc=" << p.chargeParity().parity() 
    << "\tG=" << p.gParity().parity() << "\ttheGParity=" << p.theGParity() << "\t2*I=" << p.twoIso() << "\t2*I3=" << p.twoIso3() 
    << "\tcharm=" << p.charm() << "\tstrange=" << p.strange() << std::endl;

  return o;
}
