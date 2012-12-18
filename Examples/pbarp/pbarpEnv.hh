#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/DataUtils.hh"
#include "Examples/pbarp/pbarpParser.hh"

class Particle;
class ParticleTable;
class pbarpReaction;
class IsobarDecayList;
class pbarpEventList;

class pbarpEnv {

public:
  static pbarpEnv* instance();
  ~pbarpEnv();

  void setup(pbarpParser& thePbarpParser);
  ParticleTable* particleTable() {return _particleTable;}

  const int lmax() const {return _lmax;}
  const int noFinalStateParticles() {return _noFinalStateParticles;}
  std::vector<Particle*> finalStateParticles() {return _finalStateParticles;}
  std::vector<std::pair<Particle*, Particle*> > producedParticlePairs() {return _producedParticlePairs;}
  boost::shared_ptr<IsobarDecayList> decayList() {return _decList;}
  boost::shared_ptr<IsobarDecayList> prodDecayList() {return _prodDecList;}
  boost::shared_ptr<pbarpReaction> reaction() {return _pbarpReaction;} 

protected:  

  pbarpEnv();
  static pbarpEnv* _instance;
  bool _alreadySetUp;
  int _lmax;
  int _noFinalStateParticles;

  ParticleTable* _particleTable;
  std::vector<Particle*> _finalStateParticles;
  std::vector<std::pair<Particle*, Particle*> > _producedParticlePairs;
  boost::shared_ptr<IsobarDecayList> _decList;
  boost::shared_ptr<IsobarDecayList> _prodDecList;
  boost::shared_ptr<pbarpReaction> _pbarpReaction;
};
