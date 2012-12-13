#include "Particle/Particle.hh"

#include "ErrLogger/ErrLogger.hh"

#include <cstdlib>
#include <iostream>
#include <string>
#include <cstdlib>

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

int Particle::iso()
{
  if (0 != pdata)
    return pdata->iso;
  else {
    Alert << "Error::Particle: accessing uninitialized data" << endmsg;
    exit(1);
  }
}

int Particle::iso3()
{
  if (0 != pdata)
    return pdata->isoThree;
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
      << "\tG=" << gParity().parity() << "\tI=" << iso() << "\tI3=" << iso3() 
      << "\tcharm=" << charm() << "\tstrange=" << strange() << std::endl;
  
  return;
}

Particle::Particle()
{
}


std::ostream &operator<<(std::ostream &o, Particle &p)
{
  o << p.name() << "\tmass=" << p.massErr() << "\twidth=" << p.widthErr() << "\t3*q=" << p.charge()
    << "\t2*J=" << p.twoJ() << "\tP=" << p.theParity() << "\tp=" << p.parity().parity() 
    << "\tC=" << p.theCParity()  << "\tc=" << p.chargeParity().parity() 
    << "\tG=" << p.gParity().parity() << "\tI=" << p.iso() << "\tI3=" << p.iso3() 
    << "\tcharm=" << p.charm() << "\tstrange=" << p.strange() << std::endl;

  return o;
}
