#pragma once

#include "PwaUtils/ParserBase.hh"
// Boost headers go here



class pbarpParser : public ParserBase 
{

  public:

  pbarpParser(int argc,char **argv);
  virtual ~pbarpParser(){;}

  const unsigned int getLMax() const { return _lMax; }
   const std::string pdgTableFile() const {return _pdgTableFile;} 
  const std::vector<std::string>& finalStateParticles() const { return _finalStateParticles; }
  const std::vector<std::string>& productionSystem() const { return _productionSystem; }
  const std::vector<std::string>& decaySystem() const { return _decaySystem; }
  const std::vector<std::string>& replaceSuffixNames() const { return _replaceParSuffix; }
  const std::vector<std::string>& histMassNames() const { return _histMass; }
  const std::vector<std::string>& histAngleNames() const { return _histAngles;}
  const std::vector<std::string>& decayDynamics() const { return _dynamics;}

protected:
  virtual bool parseCommandLine(int argc,char **argv); 
  unsigned int _lMax;
  std::string _pdgTableFile;  
  std::vector<std::string> _finalStateParticles;
  std::vector<std::string> _productionSystem;
  std::vector<std::string> _decaySystem;
  std::vector<std::string> _dynamics;
  std::vector<std::string> _replaceParSuffix;
  std::vector<std::string> _histMass;
  std::vector<std::string> _histAngles;
};



