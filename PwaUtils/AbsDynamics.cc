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
  ,_massKey("defaultMassKey")
  ,_fsParticles(fsParticles)
  ,_mother(mother)
  ,_dynKey(FunctionUtils::particleListName(fsParticles))
  ,_grandmaKey("default")
{
  if(0!=mother) _massKey=mother->name(); 
}

AbsDynamics::~AbsDynamics()
{
}

void AbsDynamics::cacheAmplitudes(){
  _cacheAmps=true;
}

void AbsDynamics::fillMasses(EvtData* theData){

  Vector4<double> mass4Vec(0.,0.,0.,0.);
  std::vector<Particle*>::iterator it;
  for (it=_fsParticles.begin(); it != _fsParticles.end(); ++it){
    mass4Vec+=theData->FourVecsString[(*it)->name()];
  }

  theData->FourVecsString[_dynKey]=mass4Vec;
}

