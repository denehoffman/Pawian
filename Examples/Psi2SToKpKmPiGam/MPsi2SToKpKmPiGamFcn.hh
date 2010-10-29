#ifndef _MPsi2SToKpKmPiGamFcn_H
#define _MPsi2SToKpKmPiGamFcn_H

#include <iostream>
#include <fstream>
//#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
//#include <cassert>
#include "Minuit2/FCNBase.h"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "PwaUtils/DataUtils.hh"
#include "Minuit2/MnUserParameters.h"

class AbsPsi2SToKpKmPiGamLh;


namespace ROOT {

   namespace Minuit2 {
class MPsi2SToKpKmPiGamFcn : public FCNBase {

public:

  MPsi2SToKpKmPiGamFcn(boost::shared_ptr<AbsPsi2SToKpKmPiGamLh> psi2SToKpKmPiGamLh);
  virtual ~MPsi2SToKpKmPiGamFcn();

  double operator()(const std::vector<double>& par) const;

  double Up() const;


private:
  boost::shared_ptr<AbsPsi2SToKpKmPiGamLh> _psi2SToKpKmPiGamLhPtr;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif 
