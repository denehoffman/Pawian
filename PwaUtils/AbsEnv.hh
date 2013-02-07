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
#include "PwaUtils/AbsHist.hh"

class Particle;
class ParticleTable;
class AbsDecayList;
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
  boost::shared_ptr<AbsDecayList> absDecayList() {return _absDecList;}
  boost::shared_ptr<AbsDecayList> prodDecayList() {return _prodDecList;}
  const std::string outputFileNameSuffix() const {return _outputFileNameSuffix;}
  const std::string serializationFileName() const {return _serializationFileName;}
  std::vector<std::vector<std::string> >& histMassSystems() {return _histMassSystems;} 
  std::vector<boost::shared_ptr<angleHistData> >& angleHistDataVec() {return _angleHistDataVec;} 

protected:  

  AbsEnv();
  bool _alreadySetUp;
  int _noFinalStateParticles;
  ParticleTable* _particleTable;
  std::vector<Particle*> _finalStateParticles;
  std::vector<std::pair<Particle*, Particle*> > _producedParticlePairs;
  boost::shared_ptr<AbsDecayList> _absDecList;
  boost::shared_ptr<AbsDecayList> _prodDecList;
  std::string _outputFileNameSuffix;
  std::string _serializationFileName;
  std::vector<std::vector<std::string> > _histMassSystems;
  std::vector<boost::shared_ptr<angleHistData> > _angleHistDataVec;
};
