// IsobarHeliDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/AbsEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"

IsobarHeliDecay::IsobarHeliDecay(Particle* mother, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv) :
  AbsDecay(mother, daughter1, daughter2, theEnv)
{
  validJPClamlam( _motherJPCPtr, _daughter1, _daughter2, _JPClamlamDecAmps);
}

IsobarHeliDecay::IsobarHeliDecay(boost::shared_ptr<const jpcRes> motherJPCPtr, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv, std::string motherName) :
  AbsDecay(motherJPCPtr, daughter1, daughter2, theEnv, motherName)
{
  validJPClamlam( _motherJPCPtr, _daughter1, _daughter2, _JPClamlamDecAmps);
}

IsobarHeliDecay::~IsobarHeliDecay(){
}

void IsobarHeliDecay::print(std::ostream& os) const{
  os << "\nJPClamlam amplitudes for decay\t" << _name << ":\n";
  os << "suffix for fit parameter name:\t" << _fitParamSuffix << "\n";
  
  std::vector< boost::shared_ptr<const JPClamlam> >::const_iterator it;
  for (it = _JPClamlamDecAmps.begin(); it!= _JPClamlamDecAmps.end(); ++it){
    (*it)->print(os);
    os << "\n";
  }

  AbsDecay::print(os);  
  os << "\n";
}
