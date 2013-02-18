// WoDynamics class definition file. -*- C++ -*-
// Copyright 20123Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/WoDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

WoDynamics::WoDynamics(std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsDynamics(fsParticles, mother)
{
}

WoDynamics::~WoDynamics()
{
}

complex<double> WoDynamics::eval(EvtData* theData, Spin OrbMom){

  complex<double> result(1.,0.);
  return result;
}

void  WoDynamics::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
}

bool WoDynamics::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;
  return _recalculate;
}

void WoDynamics::updateFitParams(fitParams& theParamVal){
}
