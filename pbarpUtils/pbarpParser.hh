// pbarpParser class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include "PwaUtils/ParserBase.hh"
// Boost headers go here



class pbarpParser : public ParserBase 
{

  public:

  pbarpParser(int argc,char **argv);
  virtual ~pbarpParser(){;}

  const unsigned int getLMax() const { return _lMax; }
  const float getpbarMomentum() const { return _pbarMomentum; } 
  //  const std::vector<std::string>& productionSystem() const { return _productionSystem; }
  const std::vector<std::string>& histMassNames() const { return _histMass; }
  const std::vector<std::string>& histAngleNames() const { return _histAngles;}
  const std::vector<std::string>& spinDensityNames() const { return _spinDensity;}
  const std::vector<std::string>& lmaxParticle() const{ return _lmaxParticle; }

protected:
  virtual bool parseCommandLine(int argc,char **argv); 
  unsigned int _lMax;
  float _pbarMomentum;
  //  std::vector<std::string> _productionSystem;
  std::vector<std::string> _histMass;
  std::vector<std::string> _histAngles;
  std::vector<std::string> _spinDensity;
  std::vector<std::string> _lmaxParticle;
};



