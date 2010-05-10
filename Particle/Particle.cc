#include "Particle/Particle.hh"

#include <cstdlib>
#include <iostream>
#include <string>

Particle::Particle(ParticleData& data)
{
  pdata = &data;
}

Particle::Particle(Particle& other)
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
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}

const std::string& Particle::texName()
{
  if (0 != pdata)
    return pdata->texName;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}

ParticleType Particle::type()
{
  if (0 != pdata)
    return pdata->type;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}

int Particle::charge()
{
  if (0 != pdata)
    return pdata->charge;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
  return 0;
}

const ErrValue& Particle::massErr() const
{
  if (0 != pdata)
    return pdata->mass;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}

const ErrValue& Particle::widthErr() const
{
  if (0 != pdata)
    return pdata->width;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}

double Particle::mass()
{
  if (0 != pdata)
    return (pdata->mass).mean();
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
  return 0.;
}

double Particle::width()
{
  if (0 != pdata)
    return (pdata->width).mean();
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
  return 0.;
}

int Particle::twoJ()
{
  if (0 != pdata)
    return pdata->twoJ;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
  return 0;

}

double Particle::J()
{
  if (0 != pdata)
    return double(pdata->twoJ)/2.;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
  return 0.;
}

Parity& Particle::parity()
{
  if (0 != pdata)
    return pdata->parity;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}

Parity& Particle::chargeParity()
{
  if (0 != pdata)
    return pdata->chargeParity;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}

Parity& Particle::gParity()
{
  if (0 != pdata)
    return pdata->gParity;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}

int Particle::iso()
{
  if (0 != pdata)
    return pdata->iso;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}

int Particle::iso3()
{
  if (0 != pdata)
    return pdata->isoThree;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}

int Particle::strange()
{
  if (0 != pdata)
    return pdata->strange;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}

int Particle::charm()
{
  if (0 != pdata)
    return pdata->charm;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}

DynFunctionType Particle::dynFctType()
{
  if (0 != pdata)
    return pdata->dynamicFunction;
  else {
    std::cerr << "Error::Particle: accessing uninitialized data" << std::endl;
    exit(1);
  }
}


ParticleData* Particle::data()
{
  return pdata;
}

void Particle::print(std::ostream& out)
{
  out << name() << "\tmass=" << massErr() << "\twidth=" << widthErr() << "\t3*q=" << charge()
      << "\t2*J=" << twoJ() << "\tP=" << parity().parity() << "\tC=" << chargeParity().parity() 
      << "\tG=" << gParity().parity() << "\tI=" << iso() << "\tI3=" << iso3() 
      << "\tcharm=" << charm() << "\tstrange=" << strange() << std::endl;
  
  return;
}

Particle::Particle()
{
}

