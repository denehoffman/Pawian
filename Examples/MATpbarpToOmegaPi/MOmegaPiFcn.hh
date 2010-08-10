#ifndef _MOmegaPiFcn_H
#define _MOmegaPiFcn_H

#include <iostream>
#include <fstream>
//#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
//#include <cassert>
#include "Minuit2/FCNBase.h"
#include "Examples/MATpbarpToOmegaPi/OmegaPiData.hh"
#include "Examples/MATpbarpToOmegaPi/serDataUtils.hh"
#include "Examples/MATpbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "Minuit2/MnUserParameters.h"

class OmegaPiLh;

namespace ROOT {

   namespace Minuit2 {
class MOmegaPiFcn : public FCNBase {

public:

  MOmegaPiFcn(boost::shared_ptr<OmegaPiLh> omegaPiLh);
  virtual ~MOmegaPiFcn();

  double operator()(const std::vector<double>& par) const;

  double Up() const;

  void setMnUsrParams(MnUserParameters& upar);
  void setFitParamVal(OmegaPiData::fitParamVal& theParamVal, const std::vector<double>& par) const;

private:
  boost::shared_ptr<OmegaPiLh> _omegaPiLhPtr;
  boost::shared_ptr<const pbarpToOmegaPi0States> _barpToOmegaPi0States;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif /* _etacToapi0fcn_H */
