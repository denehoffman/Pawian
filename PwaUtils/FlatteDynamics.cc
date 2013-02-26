// FlatteDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/FlatteDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/FlatteShape.hh"
#include "Utils/FunctionUtils.hh"

FlatteDynamics::FlatteDynamics(std::string& key, std::vector<Particle*>& fsParticles, Particle* mother, std::pair<Particle*, Particle*>& decPair1stChannel, std::pair<Particle*, Particle*>& decPair2ndChannel) :
  AbsDynamics(key, fsParticles, mother)
  ,_key(key)
  ,_g11Key(key+decPair1stChannel.first->name()+decPair1stChannel.second->name())
  ,_g22Key(key+decPair2ndChannel.first->name()+decPair2ndChannel.second->name())
{
  _flattePtr=boost::shared_ptr<FlatteShape>(new FlatteShape(decPair1stChannel, decPair2ndChannel));
}

FlatteDynamics::~FlatteDynamics()
{
}

complex<double> FlatteDynamics::eval(EvtData* theData, Spin OrbMom){
  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate){
    return _cachedMap[evtNo];
  }

  complex<double> result=_flattePtr->calculate(theData->FourVecsString[_dynKey].M(), _currentMass, _currentg11, _currentg22);
  
  if ( _cacheAmps){
#ifdef _OPENMP
#pragma omp critical (flatteDynCache)
    {
#endif
      _cachedMap[evtNo]=result;
#ifdef _OPENMP
    }
#endif
  }  

  return result;
}

void  FlatteDynamics::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
    fitVal.Masses[_key]=_mother->mass();
    fitErr.Masses[_key]=0.03;
    fitVal.gFactors[_g11Key]=1.;
    fitErr.gFactors[_g11Key]=1.;
    fitVal.gFactors[_g22Key]=1.;
    fitErr.gFactors[_g22Key]=1.;
}

bool FlatteDynamics::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;

  double mass=theParamVal.Masses[_key];
  if ( fabs(mass-_currentMass) > 1.e-10){
    _recalculate=true;
  }

  double g11=theParamVal.gFactors[_g11Key];
  if ( fabs(g11-_currentg11) > 1.e-10){
    _recalculate=true;
  }

  double g22=theParamVal.gFactors[_g22Key];
  if ( fabs(g22-_currentg22) > 1.e-10){
    _recalculate=true;
  }
  return _recalculate;
}

void FlatteDynamics::updateFitParams(fitParams& theParamVal){
  _currentMass=theParamVal.Masses[_key];
  _currentg11=theParamVal.gFactors[_g11Key];
  _currentg22=theParamVal.gFactors[_g22Key];
}
