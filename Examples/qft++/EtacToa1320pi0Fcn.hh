#ifndef _etacToa1320pi0fcn_H
#define _etacToa1320pi0fcn_H

#include <iostream>
#include <fstream>
//#include <string>
#include <vector>

//#include <cassert>
#include "Minuit2/FCNBase.h"


class EtacToa1320pi0fit;


namespace ROOT {

   namespace Minuit2 {
class EtacToa1320pi0Fcn : public FCNBase {

public:

  EtacToa1320pi0Fcn(EtacToa1320pi0fit* etacToa1320pi0fit);
  virtual ~EtacToa1320pi0Fcn();

  double operator()(const std::vector<double>& par) const;

  double Up() const;

private:
  EtacToa1320pi0fit* _etacToa1320pi0fit;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif /* _etacToa1320pi0fcn_H */
