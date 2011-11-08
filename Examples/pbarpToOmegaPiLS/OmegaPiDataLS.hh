#ifndef _OmegaPiDataLS_H
#define _OmegaPiDataLS_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"

namespace OmegaPiDataLS {

  struct fitParamVal
  {
    std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess > lsParam;
    pair<double, double>  omegaMass;

//     std::map< boost::shared_ptr<const JPCLSls>, boost::shared_ptr< pair<double, double> >, pawian::Collection::SharedPtrLess > singletLsParam;
//     std::map< boost::shared_ptr<const JPCLSls>, boost::shared_ptr< pair<double, double> >, pawian::Collection::SharedPtrLess > triplet0LsParam;
//     std::map< boost::shared_ptr<const JPCLSls>, boost::shared_ptr< pair<double, double> >, pawian::Collection::SharedPtrLess > tripletp1LsParam;
//     std::map< boost::shared_ptr<const JPCLSls>, boost::shared_ptr< pair<double, double> >, pawian::Collection::SharedPtrLess > tripletm1LsParam;

// //     void setMultpletParams(){
// //       std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
// //       std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itMultiplet;

// //       for ( it=lsParam.begin(); it!=lsParam.end(); ++it){
// // 	for ( itMultiplet=singletLsParam.begin(); itMultiplet!=singletLsParam.end(); ++itMultiplet){
// // 	}
// //       }
//     };
  };
  
  struct OmPiEvtDataLS
  {
    Vector4<float> cm_4Vec;
    Vector4<float> omegaHeliCm4Vec;
    Vector4<float> omegaHeliCm4Vec2;
    Vector4<float> pi0RecHeliCm4Vec;
    Vector4<float> pi0DecHeliCm4Vec;
    Vector4<float> pi0HeliOmega4Vec;
    Vector4<float> pi0HeliOmega4Vec2;

    Vector4<float> pipHeliOmega4Vec;  //for omega -> 3pi only
    Vector4<float> pimHeliOmega4Vec;  //for omega -> 3pi only
    Vector4<float> pippimCrossHeliOmega4Vec; //for omega -> 3pi only
    
    float cosPi0HeliOmega4Vec;
    float lambda; //for omega -> 3pi only
    float eventWeight;

    map<Spin,map<Spin,map<Spin,complex<double> > > > Dfp; //Wigner D functions for omega pi0 production
    map<Spin,map<Spin,map<Spin,complex<double> > > > Dfd; //Wigner D functions for omega decay to pi0 gamma
    map<Spin,map<Spin,map<Spin,complex<double> > > > Dfp2; //Wigner D functions for omega pi0 production
    map<Spin,map<Spin,map<Spin,complex<double> > > > Dfd2; //Wigner D functions for omega decay to pi0 gamma
  };
  
}
#endif /* _etacToapi0Data_H */
