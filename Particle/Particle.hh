#ifndef PARTICLE_HH
#define PARTICLE_HH

#include <string>
#include <iostream>
#include "Parity.hh"

enum ParticleType {parton, lepton, boson, meson, baryon, nucleus, special};

struct ParticleData {
  ParticleType type;
  std::string name;
  std::string texName;
  int charge;
  double mass;
  double width;
  int twoJ;
  int iso;
  int isoThree;
  int strange;
  int charm;
  Parity parity;
  Parity chargeParity;
  Parity gParity;
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
  Particle(ParticleData& data);
  Particle(Particle& other);
  ~Particle();

  const std::string& name();
  const std::string& texName();
  ParticleType type();
  int charge();
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

  ParticleData* data();

  void print(std::ostream& out);

private:
  Particle();

  ParticleData* pdata;

};

#endif
