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


private:
  boost::shared_ptr<AbsOmegaPiLh> _omegaPiLhPtr;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif /* _etacToapi0fcn_H */
