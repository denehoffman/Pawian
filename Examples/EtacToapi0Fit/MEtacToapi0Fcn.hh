#ifndef _metacToapi0fcn_H
#define _metacToapi0fcn_H

#include <iostream>
#include <fstream>
//#include <string>
#include <vector>

//#include <cassert>
#include "Minuit2/FCNBase.h"


class EtacToapi0Lh;


namespace ROOT {

   namespace Minuit2 {
class MEtacToapi0Fcn : public FCNBase {

public:

  MEtacToapi0Fcn(EtacToapi0Lh* etacToapi0Lh);
  virtual ~MEtacToapi0Fcn();

  double operator()(const std::vector<double>& par) const;

  double Up() const;

private:
  EtacToapi0Lh* _etacToapi0Lh;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif /* _etacToapi0fcn_H */
