// BreitWignerDynamics class definition file. -*- C++ -*-
// Copyright 20123Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/BreitWignerDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

BreitWignerDynamics::BreitWignerDynamics(std::string& massKey, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsDynamics(massKey, fsParticles, mother)
  ,_massKey(massKey)
{
}

BreitWignerDynamics::~BreitWignerDynamics()
{
}

complex<double> BreitWignerDynamics::eval(EvtData* theData, Spin OrbMom){
  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate){
    return _cachedMap[evtNo];
  }

  Vector4<double> mass4Vec(0.,0.,0.,0.);
  std::vector<Particle*>::iterator it;
  for (it=_fsParticles.begin(); it != _fsParticles.end(); ++it){
    mass4Vec+=theData->FourVecsString[(*it)->name()];
    }

  complex<double> result=BreitWigner(mass4Vec, _currentMass, _currentWidth);
  
  if ( _cacheAmps){
#ifdef _OPENMP
#pragma omp critical (bwDynCache)
    {
#endif
      _cachedMap[evtNo]=result;
#ifdef _OPENMP
    }
#endif
  }  

  return result;
}

void  BreitWignerDynamics::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
    fitVal.Masses[_massKey]=_mother->mass();
    fitErr.Masses[_massKey]=0.03;
    fitVal.Widths[_massKey]=_mother->width();
    fitErr.Widths[_massKey]=0.2*_mother->width();
}

bool BreitWignerDynamics::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;

  double mass=theParamVal.Masses[_massKey];
  if ( fabs(mass-_currentMass) > 1.e-10){
    _recalculate=true;
  }
  double width=theParamVal.Widths[_massKey];
  if ( fabs(width-_currentWidth) > 1.e-10){
    _recalculate=true;
  }

  return _recalculate;
}

void BreitWignerDynamics::updateFitParams(fitParams& theParamVal){
  _currentMass=theParamVal.Masses[_massKey];
  _currentWidth=theParamVal.Widths[_massKey];
}
