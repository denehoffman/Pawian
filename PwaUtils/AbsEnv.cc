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

  //final state particles
  const std::vector<std::string> finalStateParticleStr=theParser->finalStateParticles();
  
  std::vector<std::string>::const_iterator itStr;
  for ( itStr = finalStateParticleStr.begin(); itStr != finalStateParticleStr.end(); ++itStr){
    Particle* currentParticle = _particleTable->particle(*itStr);
    _finalStateParticles.push_back(currentParticle);
  }

  _noFinalStateParticles= (int) _finalStateParticles.size();

  //decays

  std::vector<std::string> decaySystem= theParser->decaySystem();
  for ( itStr = decaySystem.begin(); itStr != decaySystem.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string motherStr;
    stringStr >> motherStr;
    Particle* motherParticle = _particleTable->particle(motherStr);
    if( 0==motherParticle){
      Alert << "mother particle\t" << motherStr << "\tdoes not exist in pdtTable" << endmsg;
      exit(1);
    }
    std::string daughter1Str;
    stringStr >> daughter1Str;
    Particle* daughter1Particle = _particleTable->particle(daughter1Str);
    if( 0==daughter1Particle){
      Alert << "first daughter particle\t" << daughter1Str << "\tdoes not exist in pdtTable" << endmsg;
      exit(1);
    }

    std::string daughter2Str;
    stringStr >> daughter2Str;
    Particle* daughter2Particle = _particleTable->particle(daughter2Str);
    if( 0==daughter2Particle){
      Alert << "second daughter particle\t" << daughter2Str << "\tdoes not exist in pdtTable" << endmsg;
      exit(1);
    }

    boost::shared_ptr<IsobarDecay> tmpDec(new IsobarDecay(motherParticle, daughter1Particle, daughter2Particle));
 
    _decList->addDecay(tmpDec);
  }

  //add dynamics

  std::vector<boost::shared_ptr<IsobarDecay> > isoDecList= _decList->getList();
 
  std::vector<std::string> decDynVec = theParser->decayDynamics();
  for ( itStr = decDynVec.begin(); itStr != decDynVec.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string particleStr;
    stringStr >> particleStr;

    std::string dynStr;
    stringStr >> dynStr;

    std::vector<boost::shared_ptr<IsobarDecay> >::iterator itDyn;
    for (itDyn=isoDecList.begin(); itDyn!=isoDecList.end(); ++itDyn){
      std::string theDecName=(*itDyn)->name();
      std::string toFind=particleStr+"To";
      size_t found;
      found=theDecName.find(toFind);
      if (found!=string::npos){
	(*itDyn)->enableDynamics(dynStr);
      }
    }
  }

  //produced particle pairs
  std::vector<std::string> productionSystem = theParser->productionSystem();

  for ( itStr = productionSystem.begin(); itStr != productionSystem.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);
    std::string firstParticleStr;
    stringStr >> firstParticleStr;
    std::cout << "first particle:\t" << firstParticleStr << std::endl;

    std::string secondParticleStr;
    stringStr >> secondParticleStr;
    std::cout << "second particle:\t" << secondParticleStr << std::endl;

    Particle* firstParticle = _particleTable->particle(firstParticleStr);
    if( 0==firstParticle){
      Alert << "particle\t" << firstParticleStr << "\tdoes not exist in pdtTable" << endmsg;
      exit(1);
    }

    Particle* secondParticle = _particleTable->particle(secondParticleStr);
    if( 0==secondParticle){
      Alert << "particle\t" << secondParticleStr << "\tdoes not exist in pdtTable" << endmsg;
      exit(1);
    }

    _producedParticlePairs.push_back(make_pair(firstParticle, secondParticle));

  }


}

