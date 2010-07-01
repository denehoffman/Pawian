#ifndef _OmegaPiData_H
#define _OmegaPiData_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"

namespace OmegaPiData {

  struct fitParamVal
  {
    std::map< boost::shared_ptr<const JPCLS>, complex<double>, pawian::Collection::SharedPtrLess > omegaProd;
    std::map< boost::shared_ptr<const jpcRes>, complex<double>, pawian::Collection::SharedPtrLess > pbarpSinglet;
    std::map< boost::shared_ptr<const jpcRes>, complex<double>, pawian::Collection::SharedPtrLess > pbarpTriplet0;
    std::map< boost::shared_ptr<const jpcRes>, complex<double>, pawian::Collection::SharedPtrLess > pbarpTriplet1;
  };
  
  struct OmPiEvtData
  {
    Vector4<float> cm_4Vec;
    Vector4<float> omegaHeliCm4Vec;
    Vector4<float> pi0RecHeliCm4Vec;
    Vector4<float> pi0HeliOmega4Vec;
    //   Tensor<complex<double> > spin2_0_amp;
    //   Tensor<complex<double> > spin2_1_amp;
    //   Tensor<complex<double> > spin1_0_amp;
    //   Tensor<complex<double> > spin1_1_amp;
    
  };
  
}
#endif /* _etacToapi0Data_H */
