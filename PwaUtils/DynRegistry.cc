// DynRegistry class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/DynRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDynamics.hh"
#include "PwaUtils/BreitWignerDynamics.hh"
#include "PwaUtils/WoDynamics.hh"
#include "ErrLogger/ErrLogger.hh"


DynRegistry* DynRegistry::_instance=0;

DynRegistry* DynRegistry::instance()
{
  if (0==_instance) _instance = new DynRegistry();
  return _instance;
}

DynRegistry::DynRegistry()
{
}

DynRegistry::~DynRegistry()
{
}

boost::shared_ptr<AbsDynamics> DynRegistry::getDynamics(boost::shared_ptr<AbsDecay> theDec){

  std::string theName=theDec->massParKey();
  boost::shared_ptr<AbsDynamics> result;

  std::map<std::string, boost::shared_ptr<AbsDynamics> >::iterator it = _dynMap.find(theName);
  if (it !=_dynMap.end()){
    result=it->second;
    Info << "Dynamics\t" << theName << "\tfound!!!" << endmsg;
  }
  else{
    std::vector<Particle*> fsParticles=theDec->finalStateParticles();

    if(theDec->withDynamics()){
      if(!theDec->hasMother()){
	Alert << "no mother resonance; can not add dynamis" << endmsg;
	exit(1);
      }
      result= boost::shared_ptr<AbsDynamics>(new BreitWignerDynamics(theName, fsParticles, theDec->motherPart()));
    }

    else{
      std::string woDynName="woDynamica";
      result= boost::shared_ptr<AbsDynamics>(new WoDynamics(theName, fsParticles, theDec->motherPart())); 
    }
  }
  return result;
}

