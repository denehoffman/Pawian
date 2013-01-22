// IsobarLSDecay class definition file. -*- C++ -*-
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

class IsobarLSDecay : public AbsDecay{

public:
  IsobarLSDecay(Particle* mother, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv);
  IsobarLSDecay(boost::shared_ptr<const jpcRes> motherJPCPtr, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv, std::string motherName="pbarp");
  virtual ~IsobarLSDecay();
  //  virtual IsobarLSDecay* clone_() const = 0;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSAmps(){ return _JPCLSDecAmps;}
  virtual void print(std::ostream& os) const;
  virtual std::string type() {return "IsobarLSDecay";} 

protected:
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSDecAmps;

};
