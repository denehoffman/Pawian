// DynRegistry class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/DynRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDynamics.hh"
#include "PwaUtils/BreitWignerDynamics.hh"
#include "PwaUtils/FlatteDynamics.hh"
#include "PwaUtils/KPiSWaveDynamics.hh"
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

  std::string theName=theDec->name();
  std::string dynType=theDec->dynType();

  boost::shared_ptr<AbsDynamics> result;

  std::map<std::string, boost::shared_ptr<AbsDynamics> >::iterator it = _dynMap.find(theName);
  if (it !=_dynMap.end()){
    result=it->second;
    Info << "Dynamics\t" << theName << "\tfound!!!" << endmsg;
  }
  else{
    std::vector<Particle*> fsParticles=theDec->finalStateParticles();

    if(theDec->hasMother()){

      if(theDec->dynType()=="BreitWigner") 
    	result= boost::shared_ptr<AbsDynamics>(new BreitWignerDynamics(theName, fsParticles, theDec->motherPart()));
      else if(theDec->dynType()=="Flatte")
    	result= boost::shared_ptr<AbsDynamics>(new FlatteDynamics(theName, fsParticles, theDec->motherPart(), theDec->firstDecayChannel(), theDec->secondDecayChannel()));
      else if(theDec->dynType()=="KpiSWave")
	result= boost::shared_ptr<AbsDynamics>(new KPiSWaveDynamics(theName, fsParticles, theDec->motherPart()));
      else if(theDec->dynType()=="WoDynamics") result= boost::shared_ptr<AbsDynamics>(new WoDynamics(theName, fsParticles, theDec->motherPart()));
      else{
    	Alert << "Dyn type:\t" << theDec->dynType() << "\tdoes not exist" << endmsg;
    	exit(1);
      }
    }
    else{ // has no mother
      if(theDec->dynType()=="WoDynamics") result= boost::shared_ptr<AbsDynamics>(new WoDynamics(theName, fsParticles, theDec->motherPart()));
      else{
	Alert << "no mother resonance; can not add dynamis" 
	      << "\nDyn type:\t" << theDec->dynType() 
	      << endmsg;
	exit(1);
      }
    }
    result->setMassKey(theDec->massParKey());
    _dynMap[theName]=result;
    _dynVec.push_back(result); 
  }

  return result;
}


