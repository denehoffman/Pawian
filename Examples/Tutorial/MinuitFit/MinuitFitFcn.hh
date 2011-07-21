#ifndef _MinuitFitFcn_H
#define _MinuitFitFcn_H

#include <iostream>
#include <fstream>
//#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
//#include <cassert>
#include "Minuit2/FCNBase.h"

class MinuitFit;

namespace ROOT {

   namespace Minuit2 {
class MinuitFitFcn : public FCNBase {

public:

  MinuitFitFcn(boost::shared_ptr<MinuitFit>);
  virtual ~MinuitFitFcn();

  double operator()(const std::vector<double>& par) const;

  double Up() const;

private:
  boost::shared_ptr<MinuitFit> _minFitPtr;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif 
