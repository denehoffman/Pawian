// XdecAmpRegistry class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/OmegaTo3PiLSDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/LSDecAmps.hh"
#include "PwaUtils/LSOmegaTo3PiDecAmps.hh"
#include "ErrLogger/ErrLogger.hh"


XdecAmpRegistry* XdecAmpRegistry::_instance=0;

XdecAmpRegistry* XdecAmpRegistry::instance()
{
  if (0==_instance) _instance = new XdecAmpRegistry();
  return _instance;
}

XdecAmpRegistry::XdecAmpRegistry()
{
}

XdecAmpRegistry::~XdecAmpRegistry()
{
}

boost::shared_ptr<AbsXdecAmp> XdecAmpRegistry::getXdecAmp(boost::shared_ptr<AbsDecay> theAbsXDec){

  boost::shared_ptr<AbsXdecAmp> result;

  std::string theName=theAbsXDec->name();
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> >::iterator it = _xDecAmpMap.find(theName);
  if (it !=_xDecAmpMap.end()) result=it->second;
  else{
    if(theAbsXDec->type()=="IsobarLSDecay"){
      boost::shared_ptr<IsobarLSDecay> decLS =  boost::dynamic_pointer_cast<IsobarLSDecay>(theAbsXDec);
      result=boost::shared_ptr<AbsXdecAmp>(new LSDecAmps(decLS));
    }
    else if(theAbsXDec->type()=="OmegaTo3PiLSDecay"){
      boost::shared_ptr<OmegaTo3PiLSDecay> decOmega =  boost::dynamic_pointer_cast<OmegaTo3PiLSDecay>(theAbsXDec);
      result=boost::shared_ptr<AbsXdecAmp>(new LSOmegaTo3PiDecAmps(decOmega));
    }
    else{
      Alert << "can nor create XdecAmp object for theAbsXDec->name():\t" << theAbsXDec->name() << endmsg;
      exit(1); 
    }

    _xDecAmpMap[result->name()]=result;
  }
  return result;
}

