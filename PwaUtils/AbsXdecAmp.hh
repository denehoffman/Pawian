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

class AbsXdecAmp {

public:
  AbsXdecAmp(const std::string& name, const std::vector<std::string>& hypVec, Spin spinX=0, int parity=1);
  AbsXdecAmp(const std::string& name);
  virtual ~AbsXdecAmp();
  
  virtual complex<double> XdecAmp(Spin lamX, EvtData* theData)=0;
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr)=0;
  virtual void print(std::ostream& os) const=0;
  const std::string name() const {return _name;}
  const Spin  spinX() const {return _J_X;}
  const int  parity() const {return _parity;}
  virtual bool checkRecalculation(fitParams& theParamVal);
  virtual void cacheAmplitudes(){_cacheAmps=true;}
  virtual void updateFitParams(fitParams& theParamVal){return;}
protected:
 

  const std::string _name;
  const std::vector<std::string> _hypVec;
  Spin _J_X;
  int _parity; 
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > _currentParamMags;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > _currentParamPhis;
  double _currentXMass;
  double _currentXWidth;
  double _currentgFactorPhiPhi;
  double _currentgFactorOmegaPhi;
  bool _cacheAmps;
  bool _recalculate;
  std::map<int, std::map<Spin, complex<double> > > _cachedAmpMap;
};
