// pbarpEnv class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsEnv.hh"
#include "PwaUtils/DataUtils.hh"
#include "pbarpUtils/pbarpHist.hh"

class Particle;
class ParticleTable;
class pbarpReaction;
class IsobarDecayList;
class pbarpEventList;
class pbarpParser;
//class pbarpHist;

class pbarpEnv : public AbsEnv{

public:
  static pbarpEnv* instance();
  ~pbarpEnv();

  void setup(pbarpParser* thePbarpParser);

  const float pbarMomentum() const {return _pbarMomentum;}
  const int lmax() const {return _lmax;}

  boost::shared_ptr<pbarpReaction> reaction() {return _pbarpReaction;} 
  std::vector<std::vector<std::string> >& histMassSystems() {return _histMassSystems;} 
  std::vector<boost::shared_ptr<angleHistData> >& angleHistDataVec() {return _angleHistDataVec;} 
  std::vector<std::string>& spinDensityNames(){ return _spinDensity;}

protected:  

  pbarpEnv();
  static pbarpEnv* _instance;
  int _lmax;
  float _pbarMomentum;

  boost::shared_ptr<pbarpReaction> _pbarpReaction;
  std::vector<std::vector<std::string> > _histMassSystems;
  std::vector<boost::shared_ptr<angleHistData> > _angleHistDataVec;
  std::vector<std::string> _spinDensity;
};
