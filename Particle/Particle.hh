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
  int iso;
  int isoThree;
  int strange;
  int charm;
  Parity parity;
  Parity chargeParity;
  Parity gParity;
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
    iso(other->iso),
    isoThree(other->isoThree),
    strange(other->strange),
    charm(other->charm),
    parity(other->parity),
    gParity(other->gParity)
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
  Parity& parity();
  Parity& chargeParity();
  Parity& gParity();
  int iso();
  int iso3();
  int strange();
  int charm();
  DynFunctionType dynFctType();
  ParticleData* data() const;
  void print(std::ostream& out);

private:
  ParticleData* pdata;

friend std::ostream &operator<<(std::ostream &o, Particle &p);
};
