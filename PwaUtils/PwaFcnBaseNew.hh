#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "Minuit2/FCNBase.h"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/FitParamsBaseNew.hh"
#include "Minuit2/MnUserParameters.h"

class AbsLhNew;

namespace ROOT {
  namespace Minuit2 {
    class PwaFcnBaseNew : public FCNBase {

    public:
      PwaFcnBaseNew(boost::shared_ptr<AbsLhNew> absLh, 
		    boost::shared_ptr<FitParamsBaseNew> fitParamsBase);
      virtual ~PwaFcnBaseNew();

      double operator()(const std::vector<double>& par) const;
      double Up() const;

    private:
      boost::shared_ptr<AbsLhNew> _absLhPtr;
      boost::shared_ptr<FitParamsBaseNew> _fitParamsBasePtr;
      mutable unsigned int _fcnCounter;
      fitParamsNew _defaultFitValParms;
      fitParamsNew _defaultFitErrParms;
    };
  }  // namespace Minuit2
}  // namespace ROOT
