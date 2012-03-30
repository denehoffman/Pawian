#ifndef _MPsi2STo2K2PiGamFcn_H
#define _MPsi2STo2K2PiGamFcn_H

#include <iostream>
#include <fstream>
//#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp> 
//#include <cassert>
#include "Minuit2/FCNBase.h"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "PwaUtils/DataUtils.hh"
#include "Minuit2/MnUserParameters.h"

class AbsPsi2STo2K2PiGamLh;


namespace ROOT {
  
  namespace Minuit2 {
    class MPsi2STo2K2PiGamFcn : public FCNBase {
      
    public:
      
      MPsi2STo2K2PiGamFcn(boost::shared_ptr<AbsPsi2STo2K2PiGamLh> psi2STo2K2PiGamLh);
      virtual ~MPsi2STo2K2PiGamFcn();
      
      double operator()(const std::vector<double>& par) const;
      
      double Up() const;

//       mutable boost::mutex mutex1; 
           
    private:
      mutable boost::shared_ptr<AbsPsi2STo2K2PiGamLh> _psi2STo2K2PiGamLhPtr;
      mutable unsigned int _fcnCounter;
    };
  }  // namespace Minuit2
  
}  // namespace ROOT


#endif 
