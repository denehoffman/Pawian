// pbarpEnv class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "pbarpUtils/pbarpEnv.hh"
#include "pbarpUtils/pbarpParser.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarLSDecay.hh"
//#include "PwaUtils/IsobarDecayList.hh"
#include "pbarpUtils/pbarpReaction.hh"
//#include "pbarpUtils/pbarpEventList.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/PdtParser.hh"

pbarpEnv* pbarpEnv::_instance=0;

pbarpEnv* pbarpEnv::instance()
{
  if (0==_instance) _instance = new pbarpEnv();
  return _instance;
}

pbarpEnv::pbarpEnv() :
  AbsEnv()
  ,_lmax(0)
  ,_pbarMomentum(0)
{
}
pbarpEnv::~pbarpEnv(){
}

void pbarpEnv::setup(pbarpParser* thePbarpParser){

  AbsEnv::setup(thePbarpParser);

  //Antiproton momentum
  _pbarMomentum = thePbarpParser->getpbarMomentum();

  //Lmax
  _lmax=thePbarpParser->getLMax();

  // individual Lmax settings
  std::vector<std::string>::iterator lmaxIt;
  std::vector<std::string> theLmaxParticles = thePbarpParser->lmaxParticle();
  for(lmaxIt = theLmaxParticles.begin(); lmaxIt != theLmaxParticles.end(); ++lmaxIt){
    std::string particle;
    short lmax;
    std::stringstream stringStr;
    stringStr << (*lmaxIt);
    stringStr >> particle >> lmax;
    _lmaxParticleData[particle] = lmax;
  }


  std::vector<std::string>::const_iterator itStr;


  //pbarp reaction
  _pbarpReaction=boost::shared_ptr<pbarpReaction>(new pbarpReaction(_producedParticlePairs, _lmax));

  //fill prodDecayList
  std::vector< boost::shared_ptr<IsobarLSDecay> > prodDecs= _pbarpReaction->productionDecays();
  std::vector< boost::shared_ptr<IsobarLSDecay> >::iterator itDec;
  for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
    _prodDecList->addDecay(*itDec);
  }

  //set suffixes
  std::vector<std::string> suffixVec = thePbarpParser->replaceSuffixNames();
  std::map<std::string, std::string> decSuffixNames;

  for ( itStr = suffixVec.begin(); itStr != suffixVec.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);
    std::string classStr;
    stringStr >> classStr;

    std::string suffixStr;
    stringStr >> suffixStr;
    decSuffixNames[classStr]=suffixStr;
  }
  
  //set suffixes for decay classes
  std::map<std::string, std::string>::iterator itMapStrStr;
  for (itMapStrStr=decSuffixNames.begin(); itMapStrStr!=decSuffixNames.end(); ++itMapStrStr){
    _absDecList->replaceSuffix(itMapStrStr->first, itMapStrStr->second);
    _prodDecList->replaceSuffix(itMapStrStr->first, itMapStrStr->second);
    //    boost::shared_ptr<IsobarDecay> theDec=_decList->decay(itMapStrStr->first);
  }

  //replace mass key
  std::vector<std::string> replMassKeyVec = thePbarpParser->replaceMassKey();
  std::map<std::string, std::string> decRepMassKeyNames;

  for ( itStr = replMassKeyVec.begin(); itStr != replMassKeyVec.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);
    std::string oldStr;
    stringStr >> oldStr;

    std::string newStr;
    stringStr >> newStr;
    decRepMassKeyNames[oldStr]=newStr;
  }

  for (itMapStrStr=decRepMassKeyNames.begin(); itMapStrStr!=decRepMassKeyNames.end(); ++itMapStrStr){
    _absDecList->replaceMassKey(itMapStrStr->first, itMapStrStr->second);
  }


  //add dynamics 
  std::vector<boost::shared_ptr<AbsDecay> > absDecList= _absDecList->getList();
  std::vector<std::string> decDynVec = thePbarpParser->decayDynamics();
  for ( itStr = decDynVec.begin(); itStr != decDynVec.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string particleStr;
    stringStr >> particleStr;

    std::string dynStr;
    stringStr >> dynStr;

    std::string tmpName;
    std::vector<std::string> additionalStringVec;
    while(stringStr >> tmpName){
      additionalStringVec.push_back(tmpName);
    }

    std::vector<boost::shared_ptr<AbsDecay> >::iterator itDec;
    for (itDec=absDecList.begin(); itDec!=absDecList.end(); ++itDec){
      std::string theDecName=(*itDec)->name();
      std::string toFind=particleStr+"To";
      size_t found;
      found = theDecName.find(toFind);
      if (found!=string::npos){
	(*itDec)->enableDynamics(dynStr, additionalStringVec);
      }
    }
  }

  // std::vector<std::string> theHistMassNames=thePbarpParser->histMassNames();
  // //fill vector histMassSystems
  // for ( itStr = theHistMassNames.begin(); itStr != theHistMassNames.end(); ++itStr){
  //   std::stringstream stringStr;
  //   stringStr << (*itStr);
    
  //   std::string tmpName;
  //   std::vector<std::string> currentStringVec;
  //   while(stringStr >> tmpName){
  //     currentStringVec.push_back(tmpName);
  //   } 
  //   _histMassSystems.push_back(currentStringVec);
  // }

  std::vector<std::string> theHistAngleNames=thePbarpParser->histAngleNames();
  //fill vector histMassSystems
  for ( itStr = theHistAngleNames.begin(); itStr != theHistAngleNames.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);
    
    std::string tmpName;
    std::vector<std::string> currentStringDecVec;
    std::vector<std::string> currentStringDecVec2;
    std::vector<std::string> currentStringMotherVec;
    bool isDecParticle=true;
    bool isFirstDecParticle=true;
    short nBodyDecay=2;
    while(stringStr >> tmpName){
      if(tmpName=="from") {
	isDecParticle=false;
	continue;
      }
      else if(tmpName=="and") {
	isFirstDecParticle=false;
	nBodyDecay=3;
	continue;
      }
      if(isDecParticle && isFirstDecParticle) currentStringDecVec.push_back(tmpName);
      else if(isDecParticle && !isFirstDecParticle) currentStringDecVec2.push_back(tmpName);
      else currentStringMotherVec.push_back(tmpName);
    }
    boost::shared_ptr<angleHistData> currentAngleHistData(new angleHistData(currentStringMotherVec, currentStringDecVec, currentStringDecVec2, nBodyDecay));
    _angleHistDataVec.push_back(currentAngleHistData);
  }

  // spin density particles
  _spinDensity = thePbarpParser->spinDensityNames();

}


