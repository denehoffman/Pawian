// AbsDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/AbsParamHandler.hh"

class Particle;

class AbsDynamics : public AbsParamHandler{

public:
  AbsDynamics(std::vector<Particle*>& fsParticles, Particle* mother );
  virtual ~AbsDynamics();

  virtual complex<double> eval(EvtData* theData, Spin OrbMom=0)=0;
  virtual void cacheAmplitudes();

protected:
  std::vector<Particle*> _fsParticles;
  Particle* _mother;
  std::map<int, complex<float> >  _cachedMap;

private:

};
