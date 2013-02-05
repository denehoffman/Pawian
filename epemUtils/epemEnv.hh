// epemEnv class definition file. -*- C++ -*-
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
class epemReaction;
class IsobarDecayList;
class epemParser;

class epemEnv : public AbsEnv{

public:
  static epemEnv* instance();
  ~epemEnv();

  void setup(epemParser* theEpEmParser);

  boost::shared_ptr<epemReaction> reaction() {return _epemReaction;}
  std::vector<std::vector<std::string> >& histMassSystems() {return _histMassSystems;} 
  std::vector<boost::shared_ptr<angleHistData> >& angleHistDataVec() {return _angleHistDataVec;} 

protected:  

  epemEnv();
  static epemEnv* _instance;

  boost::shared_ptr<epemReaction> _epemReaction;
  std::vector<std::vector<std::string> > _histMassSystems;
  std::vector<boost::shared_ptr<angleHistData> > _angleHistDataVec;
  std::vector<std::string> _spinDensity;
};
