// AbsEnv class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "PwaUtils/AbsEnv.hh"
#include "PwaUtils/IsobarDecay.hh"
#include "PwaUtils/IsobarDecayList.hh"
#include "PwaUtils/ParserBase.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/PdtParser.hh"


AbsEnv::AbsEnv() :
  _alreadySetUp(false)
  , _noFinalStateParticles(0)
  ,_decList(new IsobarDecayList())
  ,_prodDecList(new IsobarDecayList())
{
}

AbsEnv::~AbsEnv(){
}

void AbsEnv::setup(ParserBase* theParser){
  if(_alreadySetUp){
    Alert << " AbsEnv already set up" << endmsg;
    exit(1);
  }
  
  _alreadySetUp=true;

  // common options (move to base class later)
  _outputFileNameSuffix = theParser->outputFileNameSuffix();

  // pdtTable
  PdtParser pdtParser;
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
  std::string pdtFileRelPath=theParser->pdgTableFile(); 
  std::string pdtFile(theSourcePath+pdtFileRelPath);
  _particleTable = new ParticleTable;
  
  if (!pdtParser.parse(pdtFile, *_particleTable)) {
    Alert << "can not parse particle table " << pdtFile << endmsg;
    exit(1);
  }


}

