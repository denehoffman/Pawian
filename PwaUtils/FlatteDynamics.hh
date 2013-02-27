// FlatteDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsDynamics.hh"

class Flatte;

class FlatteDynamics : public AbsDynamics{

public:
  FlatteDynamics(std::string& key, std::vector<Particle*>& fsParticles, Particle* mother, std::pair<Particle*, Particle*>& decPair1stChannel, std::pair<Particle*, Particle*>& decPair2ndChannel);
  virtual ~FlatteDynamics();

  virtual complex<double> eval(EvtData* theData, Spin OrbMom=0);
  
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual bool checkRecalculation(fitParams& theParamVal);
  virtual void updateFitParams(fitParams& theParamVal);

protected:
  std::string _key;
  std::string _g11Key;
  std::string _g22Key;

  double _currentMass;
  double _currentg11;
  double _currentg22;

  std::map<int, complex<double> >  _cachedMap; 
  boost::shared_ptr<Flatte> _flattePtr;

private:

};
