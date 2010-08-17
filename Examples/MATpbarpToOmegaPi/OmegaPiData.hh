#ifndef _OmegaPiData_H
#define _OmegaPiData_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <boost/archive/text_oarchive.hpp> 
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "Examples/MATpbarpToOmegaPi/serDataUtils.hh"

//#include "Examples/MATpbarpToOmegaPi/serVec4.hh"
#include "Examples/MATpbarpToOmegaPi/serSpin.hh"

namespace OmegaPiData {

  struct fitParamVal
  {
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > omegaProdSinglet;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > omegaProdTriplet0;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > omegaProdTriplet1;
  };
  
  struct OmPiEvtData
  {

    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive & ar, const unsigned int) {
      using boost::serialization::make_nvp;

      ar & make_nvp("cm_4Vec",cm_4Vec);
      ar & make_nvp("omegaHeliCm4Vec",omegaHeliCm4Vec);
      ar & make_nvp("pi0RecHeliCm4Vec",pi0RecHeliCm4Vec);
      ar & make_nvp("pi0HeliOmega4Vec",pi0HeliOmega4Vec);
      ar & make_nvp("Dfp",Dfp);
      ar & make_nvp("Dfd",Dfd);
    }

 /*   Vector4<float> cm_4Vec;
    Vector4<float> omegaHeliCm4Vec;
    Vector4<float> pi0RecHeliCm4Vec;
    Vector4<float> pi0HeliOmega4Vec;*/
    std::vector<float> cm_4Vec;
    std::vector<float> omegaHeliCm4Vec;
    std::vector<float> pi0RecHeliCm4Vec;
    std::vector<float> pi0HeliOmega4Vec;
    map<serSpin,map<serSpin,map<serSpin,complex<double> > > > Dfp; //Wigner D functions for omega pi0 production
    map<serSpin,map<serSpin,map<serSpin,complex<double> > > > Dfd; //Wigner D functions for omega decay to pi0 gamma
  };
  
}

#endif /* _etacToapi0Data_H */
