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
  AbsXdecAmp(const std::string& name, const std::vector<std::string>& hypVec);
  virtual ~AbsXdecAmp();
  
  virtual complex<double> XdecAmp(Spin lamX, EvtDataNew* theData, fitParamsNew& theParamVal)=0;
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr)=0;
  virtual void print(std::ostream& os) const=0;
  const std::string name() const {return _name;}

protected:
  const std::string _name;
  const std::vector<std::string> _hypVec;
};
