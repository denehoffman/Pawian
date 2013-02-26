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
  AbsDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother );
  virtual ~AbsDynamics();

  virtual std::string name() {return _name;};
  virtual complex<double> eval(EvtData* theData, Spin OrbMom=0)=0;
  virtual void cacheAmplitudes();
  virtual void fillMasses(EvtData* theData);

protected:
  std::string _name;
  std::vector<Particle*> _fsParticles;
  Particle* _mother;
  std::map<int, complex<float> >  _cachedMap;
  std::string _dynKey;

private:

};
