// AbsXdecAmp class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/DynRegistry.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"


AbsXdecAmp::AbsXdecAmp(boost::shared_ptr<AbsDecay> theDec) :
  AbsParamHandler()
  , _decay(theDec)
  , _name(theDec->name())
  ,_JPCPtr(theDec->motherJPC())
  ,_absDyn(theDec->getDynamics())
  ,_key("_"+theDec->fitParSuffix())
  ,_wignerDKey(theDec->wignerDKey())
  ,_daughter1IsStable(theDec->isDaughter1Stable())
  ,_daughter2IsStable(theDec->isDaughter2Stable())
{
  initialize();
}

AbsXdecAmp::~AbsXdecAmp()
{
}

void AbsXdecAmp::initialize(){

  _absDyn = DynRegistry::instance()->getDynamics(_decay);

  if(!_daughter1IsStable){
    boost::shared_ptr<AbsDecay> decDaughter1=_decay->decDaughter1();
    _decAmpDaughter1=XdecAmpRegistry::instance()->getXdecAmp(decDaughter1);
  }
  
  if(!_daughter2IsStable){
    boost::shared_ptr<AbsDecay> decDaughter2=_decay->decDaughter2();
    _decAmpDaughter2=XdecAmpRegistry::instance()->getXdecAmp(decDaughter2);
  }
  
  _Jdaughter1=(Spin) _decay->daughter1Part()->J();
  _Jdaughter2=(Spin) _decay->daughter2Part()->J();
}

complex<double> AbsXdecAmp::daughterAmp(Spin lam1, Spin lam2, EvtData* theData, Spin lamFs){
  complex<double> result(1.,0.);
  if(!_daughter1IsStable) result *= _decAmpDaughter1->XdecAmp(lam1, theData, lamFs);
  if(!_daughter2IsStable) result *= _decAmpDaughter2->XdecAmp(lam2, theData, lamFs);
  return result;
}

void AbsXdecAmp::cacheAmplitudes(){
  _cacheAmps=true;
  _absDyn->cacheAmplitudes();
  if(!_daughter1IsStable) _decAmpDaughter1->cacheAmplitudes();
  if(!_daughter2IsStable) _decAmpDaughter2->cacheAmplitudes();
}

