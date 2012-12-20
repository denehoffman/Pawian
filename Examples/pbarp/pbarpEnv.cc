#include <getopt.h>
#include <fstream>

#include "Examples/pbarp/pbarpEnv.hh"
#include "Examples/pbarp/IsobarDecay.hh"
#include "Examples/pbarp/IsobarDecayList.hh"
#include "Examples/pbarp/pbarpReaction.hh"
#include "Examples/pbarp/pbarpEventList.hh"
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
  _alreadySetUp(false)
  ,_lmax(0)
  ,_noFinalStateParticles(0)
  ,_decList(new IsobarDecayList())
  ,_prodDecList(new IsobarDecayList())
{
}
pbarpEnv::~pbarpEnv(){
}

void pbarpEnv::setup(pbarpParser& thePbarpParser){
  if(_alreadySetUp){
    Alert << " pbarpEnv already set up" << endmsg;
    exit(1);
  }

  _alreadySetUp=true;
  // pdtTable
  PdtParser pdtParser;
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
  std::string pdtFileRelPath=thePbarpParser.pdgTableFile(); 
  std::string pdtFile(theSourcePath+pdtFileRelPath);
  _particleTable = new ParticleTable;
  
  if (!pdtParser.parse(pdtFile, *_particleTable)) {
    Alert << "can not parse particle table " << pdtFile << endmsg;
    exit(1);
  }

  //Lmax
  _lmax=thePbarpParser.getLMax();

  //final state particles
  const std::vector<std::string> finalStateParticleStr=thePbarpParser.finalStateParticles();
  
  std::vector<std::string>::const_iterator itStr;
  for ( itStr = finalStateParticleStr.begin(); itStr != finalStateParticleStr.end(); ++itStr){
    Particle* currentParticle = _particleTable->particle(*itStr);
    _finalStateParticles.push_back(currentParticle);
  }

  _noFinalStateParticles= (int) _finalStateParticles.size();
  
    //decays
  std::vector<std::string> decaySystem= thePbarpParser.decaySystem();
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

  //produced particle pairs
  std::vector<std::string> productionSystem = thePbarpParser.productionSystem();

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

  _pbarpReaction=boost::shared_ptr<pbarpReaction>(new pbarpReaction(_producedParticlePairs, _lmax));

  //fill prodDecayList
  std::vector< boost::shared_ptr<IsobarDecay> > prodDecs= _pbarpReaction->productionDecays();
  std::vector< boost::shared_ptr<IsobarDecay> >::iterator itDec;
  for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
    _prodDecList->addDecay(*itDec);
  }

  //set suffixes
  std::vector<std::string> suffixVec = thePbarpParser.replaceSuffixNames();
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

  std::vector<std::string> theHistMassNames=thePbarpParser.histMassNames();
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

  std::vector<std::string> theHistAngleNames=thePbarpParser.histAngleNames();
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


