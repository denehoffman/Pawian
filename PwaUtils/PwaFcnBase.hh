#ifndef _PwaFcnBase_H
#define _PwaFcnBase_H

#include <iostream>
#include <fstream>
//#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
//#include <cassert>
#include "Minuit2/FCNBase.h"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/DataUtils.hh"
#include "Minuit2/MnUserParameters.h"

class AbsLh;
class FitParamsBase;

namespace ROOT {

   namespace Minuit2 {
class PwaFcnBase : public FCNBase {

public:

  PwaFcnBase(boost::shared_ptr<AbsLh> absLh, boost::shared_ptr<FitParamsBase> fitParamsBase);
  virtual ~PwaFcnBase();

  double operator()(const std::vector<double>& par) const;

  double Up() const;


private:
  boost::shared_ptr<AbsLh> _absLhPtr;
  boost::shared_ptr<FitParamsBase> _fitParamsBasePtr;
  unsigned int *_fcnCounter;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif 
