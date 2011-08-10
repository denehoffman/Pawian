#ifndef _MIMinuitFcn_H
#define _MIMinuitFcn_H

#include <vector>
#include <boost/shared_ptr.hpp>
//#include <cassert>
#include "Minuit2/FCNBase.h"

class MIData;

namespace ROOT {

   namespace Minuit2 {
class MIMinuitFcn : public FCNBase {

public:

  MIMinuitFcn(boost::shared_ptr<MIData> theData);
  virtual ~MIMinuitFcn();

  double operator()(const std::vector<double>& par) const;

  double Up() const;

private:
  boost::shared_ptr<MIData> _myDataPtr;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif 
