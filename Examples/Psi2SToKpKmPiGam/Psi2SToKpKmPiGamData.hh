#ifndef _Psi2SToKpKmPiGamData_H
#define _Psi2SToKpKmPiGamData_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"

namespace Psi2SToKpKmPiGamData {

  struct fitParamVal
  {
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToChiGam;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK890K;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_0_K;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_1_K;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_2_K;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa0Pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa2Pi;
    pair<double, double> BwK890;
    pair<double, double> BwK1400_0;
    pair<double, double> BwK1400_1;
    pair<double, double> BwK1400_2;
    pair<double, double> Bwa980;
    pair<double, double> Bwa2;
    double FlatMa980;
    double FlatgKK;
    double FlatgEtaPi;
  };
  
  struct Psi2SToKpKmPiGamEvtData
  {
    Vector4<float> psi2S_4V;
    Vector4<float>  chic1_HeliPsi2S_4V;
    Vector4<float> KpKm_HeliChic1_4V;
    Vector4<float> KpPi_HeliChic1_4V;
    Vector4<float> KmPi_HeliChic1_4V;
    Vector4<float> Kp_HeliKpKm_4V;
    Vector4<float> Kp_HeliKpPi_4V;
    Vector4<float> Km_HeliKmPi_4V;
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfPsi; //Wigner D functions for Psi(2S) decay
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfChiToKpPi; //Wigner D functions for Chi_c1 -> (K+ pi) K-
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfChiToKmPi; //Wigner D functions for Chi_c1 -> (K- pi) K+
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfChiToa0Pi; //Wigner D functions for Chi_c1 ->  a(980) pi
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfChiToa2Pi; //Wigner D functions for Chi_c1 ->  a2 pi

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKp1; //Wigner D functions for K*1+ decay
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKp2; //Wigner D functions for K*2+ decay
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKm1; //Wigner D functions for K*1- decay
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKm2; //Wigner D functions for K*2- decay
    map<Spin,map<Spin,map<Spin,complex<double> > > > Dfa2; //Wigner D functions for a2 decay
  };
  
}
#endif /* _Psi2SToKpKmPiGamData_H */
