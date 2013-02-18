// BreitWignerDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsDynamics.hh"

class BreitWignerDynamics : public AbsDynamics{

public:
  BreitWignerDynamics(std::string& massKey, std::vector<Particle*>& fsParticles, Particle* mother);
  virtual ~BreitWignerDynamics();

  virtual complex<double> eval(EvtData* theData, Spin OrbMom=0);
  
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual bool checkRecalculation(fitParams& theParamVal);
  virtual void updateFitParams(fitParams& theParamVal);

protected:
  std::string _massKey;

  double _currentMass;
  double _currentWidth;
  std::map<int, complex<double> >  _cachedMap; 

private:

};
