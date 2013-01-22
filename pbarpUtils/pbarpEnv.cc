// pbarpEnv class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "pbarpUtils/pbarpEnv.hh"
#include "pbarpUtils/pbarpParser.hh"
#include "PwaUtils/IsobarDecay.hh"
#include "PwaUtils/IsobarDecayList.hh"
#include "pbarpUtils/pbarpReaction.hh"
#include "pbarpUtils/pbarpEventList.hh"
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

  std::vector<std::string>::const_iterator itStr;


  //pbarp reaction
  _pbarpReaction=boost::shared_ptr<pbarpReaction>(new pbarpReaction(_producedParticlePairs, _lmax));

  //fill prodDecayList
  std::vector< boost::shared_ptr<IsobarDecay> > prodDecs= _pbarpReaction->productionDecays();
  std::vector< boost::shared_ptr<IsobarDecay> >::iterator itDec;
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
    _decList->replaceSuffix(itMapStrStr->first, itMapStrStr->second);
    _prodDecList->replaceSuffix(itMapStrStr->first, itMapStrStr->second);
    boost::shared_ptr<IsobarDecay> theDec=_decList->decay(itMapStrStr->first);
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
    _decList->replaceMassKey(itMapStrStr->first, itMapStrStr->second);
  }



  std::vector<std::string> theHistMassNames=thePbarpParser->histMassNames();
  //fill vector histMassSystems
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

  std::vector<std::string> theHistAngleNames=thePbarpParser->histAngleNames();
  //fill vector histMassSystems
  for ( itStr = theHistAngleNames.begin(); itStr != theHistAngleNames.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);
    
    std::string tmpName;
    std::vector<std::string> currentStringDecVec;
    std::vector<std::string> currentStringMotherVec;
    bool isDecParticle=true;
    while(stringStr >> tmpName){
      if(tmpName=="from") {
	isDecParticle=false;
	continue;
      }
      if(isDecParticle) currentStringDecVec.push_back(tmpName);
      else currentStringMotherVec.push_back(tmpName);
    }
    boost::shared_ptr<angleHistData> currentAngleHistData(new angleHistData(currentStringMotherVec, currentStringDecVec)); 
    _angleHistDataVec.push_back(currentAngleHistData);
  }


}


