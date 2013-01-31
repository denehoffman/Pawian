// AbsXdecAmp class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FitParamsBase.hh"

class AbsDecay;

class AbsXdecAmp {

public:
  AbsXdecAmp(boost::shared_ptr<AbsDecay> theDec);
  virtual ~AbsXdecAmp();
  
  virtual complex<double> XdecAmp(Spin lamX, EvtData* theData, Spin lamFs)=0;
  virtual complex<double> XdecPartAmp(Spin lamX, Spin lamDec, short fixDaughterNr,
                                      EvtData* theData, Spin lamFs)=0;
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr)=0;
  virtual complex<double> daughterAmp(Spin lam1, Spin lam2, EvtData* theData, Spin lamFs);
  virtual void print(std::ostream& os) const=0;
  const std::string name() const {return _name;}
  const Spin  spinX() const {return _J_X;}
  const int  parity() const {return _parity;}
  virtual bool checkRecalculation(fitParams& theParamVal);
  virtual boost::shared_ptr<const jpcRes>& jpcPtr() {return _JPCPtr;}
  boost::shared_ptr<AbsDecay> absDec() {return _decay;}
  virtual void cacheAmplitudes(){_cacheAmps=true;}
  virtual void updateFitParams(fitParams& theParamVal){return;}
protected:
 
  boost::shared_ptr<AbsDecay> _decay; 
  const std::string _name;
  boost::shared_ptr<const jpcRes> _JPCPtr;
  const std::vector<std::string> _hypVec;
  Spin _J_X;
  int _parity;
  boost::shared_ptr<AbsXdecAmp> _decAmpDaughter1;
  boost::shared_ptr<AbsXdecAmp> _decAmpDaughter2; 
  Spin _Jdaughter1;
  Spin _Jdaughter2;
  std::string _key;
  std::string _massKey;
  const std::string _wignerDKey;
  bool _withDyn;

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > _currentParamMags;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > _currentParamPhis;
  double _currentXMass;
  double _currentXWidth;
  double _currentgFactorPhiPhi;
  double _currentgFactorOmegaPhi;
  bool _daughter1IsStable;
  bool _daughter2IsStable;
  bool _cacheAmps;
  bool _recalculate;
  std::map<int, std::map<Spin, std::map<Spin, complex<double> > > > _cachedAmpMap;

  virtual void initialize();
};
