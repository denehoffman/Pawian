// OmegaTo3PiLSDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

class Particle;
class EvtData;
class AbsEnv;

class OmegaTo3PiLSDecay : public AbsDecay{

public:
  OmegaTo3PiLSDecay(Particle* mother, Particle* daughter1, Particle* daughter2, Particle* daughter3, AbsEnv* theEnv);
  virtual ~OmegaTo3PiLSDecay();
  //  virtual OmegaTo3PiLSDecay* clone_() const = 0;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSAmps(){ return _JPCLSDecAmps;}
  virtual void fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, EvtData* evtData);
  virtual void print(std::ostream& os) const;
  Particle* daughter3Part() {return _daughter3;}
  virtual std::string type() {return "OmegaTo3PiLSDecay";} 

protected:
  Particle* _daughter3;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSDecAmps;
  std::vector<Particle*> _finalStateParticlesDaughter3;

};
