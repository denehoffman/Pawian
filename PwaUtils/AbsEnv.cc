// AbsEnv class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "PwaUtils/AbsEnv.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/OmegaTo3PiLSDecay.hh"
#include "PwaUtils/ParserBase.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/PdtParser.hh"


AbsEnv::AbsEnv() :
  _alreadySetUp(false)
  , _noFinalStateParticles(0)
  ,_absDecList(new AbsDecayList())
  ,_prodDecList(new AbsDecayList())
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

  // common options
  _outputFileNameSuffix = theParser->outputFileNameSuffix();
  _serializationFileName = theParser->serializationFile();

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

  // cloned particles
  const std::vector<std::string> cloneParticle=theParser->cloneParticle();
  std::vector<std::string>::const_iterator itcP;

  for ( itcP = cloneParticle.begin(); itcP != cloneParticle.end(); ++itcP){
     std::istringstream particles(*itcP);
     std::string particleOld;
     std::string particleNew;
     particles >> particleOld >> particleNew;

     _particleTable->clone(particleNew, particleOld);
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

    Particle* motherParticle =0;
    std::vector<Particle*> daughterParticles;

    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string tmpName;

    bool isDecParticle=false;
    bool firstArgument=true;
    std::string usedSystem("");

    while(stringStr >> tmpName){
      if(firstArgument){
	usedSystem=tmpName;
	firstArgument=false;
	continue;
      }
      if(tmpName=="To") {
        isDecParticle=true;
        continue;
      }
      if(isDecParticle){
  	daughterParticles.push_back(_particleTable->particle(tmpName));
      }
      else{
  	motherParticle = _particleTable->particle(tmpName);
      }
    }
    boost::shared_ptr<AbsDecay> tmpDec;
    if(daughterParticles.size()==2){
      if (usedSystem=="Heli") tmpDec= boost::shared_ptr<AbsDecay>(new IsobarHeliDecay(motherParticle, daughterParticles[0], daughterParticles[1], this));
      else if (usedSystem=="Cano")  tmpDec= boost::shared_ptr<AbsDecay>(new IsobarLSDecay(motherParticle, daughterParticles[0], daughterParticles[1], this));
      else {
	Alert << "used decay system\t" << usedSystem << "\tnot supported!!!\n" << endmsg;
	exit(1);
      }

    }

    else if(daughterParticles.size()==3){
      if (usedSystem=="Cano") tmpDec= boost::shared_ptr<AbsDecay>(new OmegaTo3PiLSDecay(motherParticle, daughterParticles[0], daughterParticles[1], daughterParticles[2], this));
      else {
	Alert << "used decay system\t" << usedSystem << "\tnot supported!!!\n" << endmsg;
	exit(1);
      }
    }
 
    else {
      Alert << "Decay\t" << (*itStr) << "\tnot supported!!!" ; 
    }

    //    _decList->addDecay(tmpDec);
    _absDecList->addDecay(tmpDec);
  }

  // //set suffixes
  // std::vector<std::string> suffixVec = theParser->replaceSuffixNames();
  // std::map<std::string, std::string> decSuffixNames;

  // //  std::vector<std::string>::iterator itStr;
  // for ( itStr = suffixVec.begin(); itStr != suffixVec.end(); ++itStr){
  //   std::stringstream stringStr;
  //   stringStr << (*itStr);
  //   std::string classStr;
  //   stringStr >> classStr;

  //   std::string suffixStr;
  //   stringStr >> suffixStr;
  //   decSuffixNames[classStr]=suffixStr;
  // }
  
  // //set suffixes for decay classes
  // std::map<std::string, std::string>::iterator itMapStrStr;
  // for (itMapStrStr=decSuffixNames.begin(); itMapStrStr!=decSuffixNames.end(); ++itMapStrStr){
  //   _absDecList->replaceSuffix(itMapStrStr->first, itMapStrStr->second);
  //   //    _prodDecList->replaceSuffix(itMapStrStr->first, itMapStrStr->second);
  //   //    boost::shared_ptr<IsobarDecay> theDec=_decList->decay(itMapStrStr->first);
  // }

  // //replace mass key
  // std::vector<std::string> replMassKeyVec = theParser->replaceMassKey();
  // std::map<std::string, std::string> decRepMassKeyNames;

  // for ( itStr = replMassKeyVec.begin(); itStr != replMassKeyVec.end(); ++itStr){
  //   std::stringstream stringStr;
  //   stringStr << (*itStr);
  //   std::string oldStr;
  //   stringStr >> oldStr;

  //   std::string newStr;
  //   stringStr >> newStr;
  //   decRepMassKeyNames[oldStr]=newStr;
  // }

  // for (itMapStrStr=decRepMassKeyNames.begin(); itMapStrStr!=decRepMassKeyNames.end(); ++itMapStrStr){
  //   _absDecList->replaceMassKey(itMapStrStr->first, itMapStrStr->second);
  // }

  // //add dynamics 
  // std::vector<boost::shared_ptr<AbsDecay> > absDecList= _absDecList->getList();
  // std::vector<std::string> decDynVec = theParser->decayDynamics();
  // for ( itStr = decDynVec.begin(); itStr != decDynVec.end(); ++itStr){
  //   std::stringstream stringStr;
  //   stringStr << (*itStr);

  //   std::string particleStr;
  //   stringStr >> particleStr;

  //   std::string dynStr;
  //   stringStr >> dynStr;

  //   std::string tmpName;
  //   std::vector<std::string> additionalStringVec;
  //   while(stringStr >> tmpName){
  //     additionalStringVec.push_back(tmpName);
  //   }

  //   std::vector<boost::shared_ptr<AbsDecay> >::iterator itDec;
  //   for (itDec=absDecList.begin(); itDec!=absDecList.end(); ++itDec){
  //     std::string theDecName=(*itDec)->name();
  //     std::string toFind=particleStr+"To";
  //     size_t found;
  //     found = theDecName.find(toFind);
  //     if (found!=string::npos){
  // 	(*itDec)->enableDynamics(dynStr, additionalStringVec);
  //     }
  //   }
  // }

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


  //fill vector histMassSystems
  std::vector<std::string> theHistMassNames=theParser->histMassNames();
  for ( itStr = theHistMassNames.begin(); itStr != theHistMassNames.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);
    
    std::string tmpName;
    std::vector<std::string> currentStringVec;
    while(stringStr >> tmpName){
      currentStringVec.push_back(tmpName);
    } 
    _histMassSystems.push_back(currentStringVec);
  }


}

