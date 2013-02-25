// AbsDynamics class definition file. -*- C++ -*-
// Copyright 20123Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/AbsDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/FunctionUtils.hh"

AbsDynamics::AbsDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsParamHandler()
  ,_name(name)
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

