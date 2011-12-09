#ifndef _MJpsiGamKsKlKKFcn_H
#define _MJpsiGamKsKlKKFcn_H

#include <iostream>
#include <fstream>
//#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
//#include <cassert>
#include "Minuit2/FCNBase.h"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKData.hh"
#include "PwaUtils/DataUtils.hh"
#include "Minuit2/MnUserParameters.h"

class AbsJpsiGamKsKlKKLh;


namespace ROOT {

   namespace Minuit2 {
class MJpsiGamKsKlKKFcn : public FCNBase {

public:

  MJpsiGamKsKlKKFcn(boost::shared_ptr<AbsJpsiGamKsKlKKLh> psi2SToKpKmPiGamLh);
  virtual ~MJpsiGamKsKlKKFcn();

  double operator()(const std::vector<double>& par) const;

  double Up() const;


private:
  boost::shared_ptr<AbsJpsiGamKsKlKKLh> _JpsiGamKsKlKKLhPtr;
  unsigned int *_fcnCounter;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif 
