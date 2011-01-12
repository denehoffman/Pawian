#ifndef _MOmegaPiFcn_H
#define _MOmegaPiFcn_H

#include <iostream>
#include <fstream>
//#include <string>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "Examples/pbarpToOmegaPi/minuitstartparam.hh"


//#include <cassert>
#include "Minuit2/FCNBase.h"
#include "Examples/pbarpToOmegaPi/OmegaPiData.hh"
#include "PwaUtils/DataUtils.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "Minuit2/MnUserParameters.h"
using namespace std;

class AbsOmegaPiLh;

namespace ROOT {

   namespace Minuit2 {
class MOmegaPiFcn : public FCNBase {

public:

  MOmegaPiFcn(boost::shared_ptr<AbsOmegaPiLh> absOmegaPiLh);
  virtual ~MOmegaPiFcn();

  double operator()(const std::vector<double>& par) const;

  double Up() const;

  void setMnUsrParams(MnUserParameters& upar);
  void setMnUsrParams(MnUserParameters& upar, minuitStartParam &theStartParam);
  void setMnUsrParams(MnUserParameters& upar, OmegaPiData::fitParamVal &finalFitParm);
  void setFitParamVal(OmegaPiData::fitParamVal& theParamVal, const std::vector<double>& par) const;

private:
  boost::shared_ptr<AbsOmegaPiLh> _omegaPiLhPtr;
  boost::shared_ptr<const pbarpToOmegaPi0States> _barpToOmegaPi0States;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif /* _etacToapi0fcn_H */
