#ifndef _MOmegaPiFcnLS_H
#define _MOmegaPiFcnLS_H

#include <iostream>
#include <fstream>
//#include <string>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "Examples/pbarpToOmegaPiLS/MinuitstartparamLS.hh"


//#include <cassert>
#include "Minuit2/FCNBase.h"

using namespace std;

class AbsOmegaPiLhLS;

namespace ROOT {

   namespace Minuit2 {
class MOmegaPiFcnLS : public FCNBase {

public:

  MOmegaPiFcnLS(boost::shared_ptr<AbsOmegaPiLhLS> absOmegaPiLh);
  virtual ~MOmegaPiFcnLS();

  double operator()(const std::vector<double>& par) const;

  double Up() const;


private:
  boost::shared_ptr<AbsOmegaPiLhLS> _omegaPiLhPtr;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif /* _etacToapi0fcn_H */
