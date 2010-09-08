#ifndef _OmegaPiData_H
#define _OmegaPiData_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
//#include "Examples/MATpbarpToOmegaPi/serDataUtils.hh"
#include "PwaUtils/DataUtils.hh"

//#include "Examples/MATpbarpToOmegaPi/serVec4.hh"
//#include "Examples/MATpbarpToOmegaPi/serSpin.hh"

namespace OmegaPiData {

  struct fitParamVal
  {
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > omegaProdSinglet;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > omegaProdTriplet0;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > omegaProdTriplet1;
  };
  
  struct OmPiEvtData
  {

    Vector4<float> cm_4Vec;
    Vector4<float> omegaHeliCm4Vec;
    Vector4<float> pi0RecHeliCm4Vec;
    Vector4<float> pi0HeliOmega4Vec;
    map<Spin,map<Spin,map<Spin,complex<double> > > > Dfp; //Wigner D functions for omega pi0 production
    map<Spin,map<Spin,map<Spin,complex<double> > > > Dfd; //Wigner D functions for omega decay to pi0 gamma
  };
  
}

#endif /* _etacToapi0Data_H */
