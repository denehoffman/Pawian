#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "Minuit2/FCNBase.h"
//#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "Minuit2/MnUserParameters.h"

class AbsLh;

namespace ROOT {
  namespace Minuit2 {
    class PwaFcnBase : public FCNBase {

    public:
      PwaFcnBase(boost::shared_ptr<AbsLh> absLh, 
		    boost::shared_ptr<FitParamsBase> fitParamsBase, std::string suffix="");
      virtual ~PwaFcnBase();

      double operator()(const std::vector<double>& par) const;
      double Up() const;

    private:
      boost::shared_ptr<AbsLh> _absLhPtr;
      boost::shared_ptr<FitParamsBase> _fitParamsBasePtr;
      mutable unsigned int _fcnCounter;
      fitParams _defaultFitValParms;
      fitParams _defaultFitErrParms;
      std::string _currentResFileName;
    };
  }  // namespace Minuit2
}  // namespace ROOT
