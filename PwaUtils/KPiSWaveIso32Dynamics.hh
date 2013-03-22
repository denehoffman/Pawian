// KPiSWaveIso32Dynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsDynamics.hh"

class AbsXdecAmp;
class KMatrixKPiSFocus;
class FVector;
class PVectorKPiSFocus;

class KPiSWaveIso32Dynamics : public AbsDynamics{

public:
  KPiSWaveIso32Dynamics(std::string& massKey, std::vector<Particle*>& fsParticles, Particle* mother);
  virtual ~KPiSWaveIso32Dynamics();

  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom=0);
  
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual bool checkRecalculation(fitParams& theParamVal);
  virtual void updateFitParams(fitParams& theParamVal);
  virtual void addGrandMa(boost::shared_ptr<AbsDecay> theDec);
  virtual const std::string& grandMaKey(AbsXdecAmp* grandmaAmp);

protected:
  //  std::string _massKey;
  double _currentMass;
  std::map<int, complex<double> >  _cachedMap;
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> >  _grandMaAmpMap;

  boost::shared_ptr<KMatrixKPiSFocus> _kMatr;
  std::map<std::string, boost::shared_ptr<FVector> > _fVecMap;
  std::map<std::string, boost::shared_ptr<PVectorKPiSFocus> > _pVecMap;

  //  std::map<std::string, std::map<std::string, double> > _currentbFactorMap;
  std::map<std::string, std::map<std::string, double> > _currentaProdMap; 
  std::map<std::string, std::map<std::string, double> > _currentbProdMap; 
  std::map<std::string, std::map<std::string, double> > _currentcProdMap;
  std::map<std::string, std::map<std::string, double> > _currentphaseProdMap;

  std::map<int, std::map<std::string, complex<float> > > _cachedStringMap;
  std::map<std::string, bool > _recalcMap;
private:

};
