// WoDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsDynamics.hh"


class WoDynamics : public AbsDynamics{

public:
  WoDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother);
  virtual ~WoDynamics();

  virtual complex<double> eval(EvtData* theData, Spin OrbMom=0);
  
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual bool checkRecalculation(fitParams& theParamVal);
  virtual void updateFitParams(fitParams& theParamVal);

protected:

private:

};
