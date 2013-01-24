// LSOmegaTo3PiDecAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/LSOmegaTo3PiDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/OmegaTo3PiLSDecay.hh"
//#include "PwaUtils/XdecAmpRegistry.hh"
#include "Particle/Particle.hh"

#ifdef _OPENMP
#include <omp.h>
#endif

LSOmegaTo3PiDecAmps::LSOmegaTo3PiDecAmps(boost::shared_ptr<OmegaTo3PiLSDecay> theDec) :
  LSDecAmps(theDec)
{
  _JPCLSs=theDec->JPCLSAmps();
  _factorMag=1.;
  if(_JPCLSs.size()>0) _factorMag=1./sqrt(_JPCLSs.size());
}

LSOmegaTo3PiDecAmps::~LSOmegaTo3PiDecAmps()
{
}


complex<double> LSOmegaTo3PiDecAmps::XdecPartAmp(Spin lamX, Spin lamDec, short fixDaughterNr, EvtData* theData, Spin lamFs){

  complex<double> result(0.,0.);
  std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
  for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
    if( fabs(lamX) > (*it)->J ) continue;
    double theMag=_currentParamMags[*it];
    double thePhi=_currentParamPhis[*it];
    complex<double> expi(cos(thePhi), sin(thePhi));

        complex<double> amp = theMag*expi*sqrt(2*(*it)->L+1)
           *conj( theData->WignerDsString[_wignerDKey][(*it)->J][lamX][0]);
        result+=amp;
  }

  return result;
}




complex<double> LSOmegaTo3PiDecAmps::XdecAmp(Spin lamX, EvtData* theData, Spin lamFs){
int evtNo=theData->evtNo;
  
  if ( _cacheAmps && !_recalculate){
    complex<double> result(0.,0.);
    result= _cachedAmpMap[evtNo][lamX][lamFs];
    return result;
  }

  complex<double> result(0.,0.);
  std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
  for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
    if( fabs(lamX) > (*it)->J ) continue;
    double theMag=_currentParamMags[*it];
    double thePhi=_currentParamPhis[*it];
    complex<double> expi(cos(thePhi), sin(thePhi));

    complex<double> amp = theMag*expi*sqrt(2*(*it)->L+1)
      *conj( theData->WignerDsString[_wignerDKey][(*it)->J][lamX][0]);
    
    result+=amp;
  }
  result*=sqrt( theData->DoubleString["lamOmegaDec"] );
  if(_withDyn){
    Vector4<double> mass4Vec(0.,0.,0.,0.);
    std::vector<Particle*> fsParticleVec=_decay->finalStateParticles();

    std::vector<Particle*>::iterator itPartVec;
    for (itPartVec=fsParticleVec.begin(); itPartVec!=fsParticleVec.end(); ++itPartVec){
      mass4Vec+=theData->FourVecsString[(*itPartVec)->name()];
    }
    result*=BreitWigner(mass4Vec, _currentXMass, _currentXWidth);
  }

  if ( _cacheAmps){
#ifdef _OPENMP
#pragma omp critical
    {
#endif
      _cachedAmpMap[evtNo][lamX][lamFs]=result;
#ifdef _OPENMP
    }
#endif
  }
  return result;
}


void LSOmegaTo3PiDecAmps::print(std::ostream& os) const{
  return; //dummy
}


