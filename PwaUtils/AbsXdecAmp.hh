// AbsXdecAmp class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <mutex>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/AbsParamHandler.hh"
#include "PwaUtils/AbsDynamics.hh"

class AbsDecay;
//class AbsDynamics;

typedef boost::unordered_map<int, std::map<Spin, std::map<Spin, complex<float> > > > intSpinSpinFloatUsMap;
typedef boost::unordered_map<std::string, intSpinSpinFloatUsMap> stringIntSpiSpinFloatMap; 

class AbsXdecAmp : public AbsParamHandler{

public:
  AbsXdecAmp(boost::shared_ptr<AbsDecay> theDec);
  virtual ~AbsXdecAmp();
  
  virtual complex<double> XdecAmp(Spin lamX, EvtData* theData, Spin lamFs, AbsXdecAmp* grandmaAmp=0)=0;
  virtual complex<double> XdecPartAmp(Spin lamX, Spin lamDec, short fixDaughterNr,
                                      EvtData* theData, Spin lamFs, AbsXdecAmp* grandmaAmp=0)=0;

  virtual complex<double> daughterAmp(Spin lam1, Spin lam2, EvtData* theData, Spin lamFs, AbsXdecAmp* grandmaAmp=0); 
  virtual void print(std::ostream& os) const=0;
  const std::string name() const {return _name;}
  virtual boost::shared_ptr<const jpcRes>& jpcPtr() {return _JPCPtr;}
  boost::shared_ptr<AbsDecay> absDec() {return _decay;}

  virtual void cacheAmplitudes();

protected:
 
  boost::shared_ptr<AbsDecay> _decay; 
  const std::string _name;
  boost::shared_ptr<const jpcRes> _JPCPtr;
  boost::shared_ptr<AbsDynamics> _absDyn;
  const std::vector<std::string> _hypVec;
  boost::shared_ptr<AbsXdecAmp> _decAmpDaughter1;
  boost::shared_ptr<AbsXdecAmp> _decAmpDaughter2; 
  Spin _Jdaughter1;
  Spin _Jdaughter2;
  std::string _key;
  const std::string _wignerDKey;
  std::mutex theMutex;

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > _currentParamMags;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > _currentParamPhis;

  bool _daughter1IsStable;
  bool _daughter2IsStable;
  bool _enabledlamFsDaughter1;
  bool _enabledlamFsDaughter2;

  //  intSpinSpinFloatUsMap _cachedAmpMap;
  stringIntSpiSpinFloatMap _cachedGrandmaAmpMap;
  virtual void initialize();
};
