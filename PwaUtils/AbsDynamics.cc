// AbsDynamics class definition file. -*- C++ -*-
// Copyright 20123Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/AbsDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

AbsDynamics::AbsDynamics(std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsParamHandler()
  ,_fsParticles(fsParticles)
  ,_mother(mother)
{
}

AbsDynamics::~AbsDynamics()
{
}

void AbsDynamics::cacheAmplitudes(){
  _cacheAmps=true;
}

