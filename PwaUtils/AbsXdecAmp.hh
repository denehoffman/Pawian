#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/EvtDataBaseListNew.hh"
#include "PwaUtils/FitParamsBaseNew.hh"

class AbsXdecAmp {

public:
  AbsXdecAmp(const std::string& name, const std::vector<std::string>& hypVec, Spin spinX=0, int parity=1);
  virtual ~AbsXdecAmp();
  
  virtual complex<double> XdecAmp(Spin lamX, EvtDataNew* theData, fitParamsNew& theParamVal)=0;
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr)=0;
  virtual void print(std::ostream& os) const=0;
  const std::string name() const {return _name;}
  const Spin  spinX() const {return _J_X;}
  const int  parity() const {return _parity;}
  virtual void checkRecalculation(fitParamsNew& theParamVal);
  void cacheAmplitudes(){_cacheAmps=true;}
protected:
  const std::string _name;
  const std::vector<std::string> _hypVec;
  Spin _J_X;
  int _parity; 
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > _oldParamMags;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > _oldParamPhis;
  double _oldXMass;
  double _oldXWidth;
  double _oldgFactorPhiPhi;
  double _oldgFactorOmegaPhi;
  bool _cacheAmps;
  bool _recalculate;
  std::map<int, std::map<Spin, complex<double> > > _cachedAmpMap;
};
