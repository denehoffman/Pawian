// epemEnv class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "epemUtils/epemEnv.hh"
#include "epemUtils/epemParser.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "epemUtils/epemReaction.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/PdtParser.hh"

epemEnv* epemEnv::_instance=0;

epemEnv* epemEnv::instance()
{
  if (0==_instance) _instance = new epemEnv();
  return _instance;
}

epemEnv::epemEnv() :
  AbsEnv()
{
}
epemEnv::~epemEnv(){
}

void epemEnv::setup(epemParser* theEpEmParser){

  AbsEnv::setup(theEpEmParser);

  std::vector<std::string>::const_iterator itStr;


  //epem reaction
  _epemReaction=boost::shared_ptr<epemReaction>(new epemReaction(_producedParticlePairs));

  //fill prodDecayList
  std::vector< boost::shared_ptr<IsobarLSDecay> > prodDecs= _epemReaction->productionCanoDecays();
  std::vector< boost::shared_ptr<IsobarLSDecay> >::iterator itDec;
  for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
    _prodDecList->addDecay(*itDec);
  }

  //set suffixes
  std::vector<std::string> suffixVec = theEpEmParser->replaceSuffixNames();
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
  std::vector<std::string> replMassKeyVec = theEpEmParser->replaceMassKey();
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



  // std::vector<std::string> theHistMassNames=theEpEmParser->histMassNames();
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

  std::vector<std::string> theHistAngleNames=theEpEmParser->histAngleNames();
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


}


