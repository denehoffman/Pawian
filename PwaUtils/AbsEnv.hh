// AbsEnv class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/DataUtils.hh"

class Particle;
class ParticleTable;
class IsobarDecayList;
class ParserBase;

class AbsEnv {

public:
  ~AbsEnv();

  void setup(ParserBase* theParser);

  ParticleTable* particleTable() {return _particleTable;}

  const int noFinalStateParticles() {return _noFinalStateParticles;}
  std::vector<Particle*> finalStateParticles() {return _finalStateParticles;}
  std::vector<std::pair<Particle*, Particle*> > producedParticlePairs() {return _producedParticlePairs;}
  boost::shared_ptr<IsobarDecayList> decayList() {return _decList;}
  boost::shared_ptr<IsobarDecayList> prodDecayList() {return _prodDecList;}
  const std::string outputFileNameSuffix() const {return _outputFileNameSuffix;}

protected:  

  AbsEnv();
  bool _alreadySetUp;
  int _noFinalStateParticles;
  ParticleTable* _particleTable;
  std::vector<Particle*> _finalStateParticles;
  std::vector<std::pair<Particle*, Particle*> > _producedParticlePairs;
  boost::shared_ptr<IsobarDecayList> _decList;
  boost::shared_ptr<IsobarDecayList> _prodDecList;
  //  std::vector<std::vector<std::string> > _histMassSystems;
  //  std::vector<boost::shared_ptr<angleHistData> > _angleHistDataVec;
  std::string _outputFileNameSuffix;
};
