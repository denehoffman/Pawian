#ifndef _Psi2STo2K2PiGamData_H
#define _Psi2STo2K2PiGamData_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"


namespace Psi2STo2K2PiGamData {

  struct fitParamVal
  {
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToChiGam;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2K892;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2K_2_1430;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2K_0_1430;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400ToK892pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K1400ToK892Pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_1270ToK892Pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_1270ToK_0_1430Pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_pif1710_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_kf1710_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980f980;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980f2200;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf1710;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_pif1370_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_kf1370_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_2400ToKf980;

    pair<double, double> BwK892;
    pair<double, double> BwK_2_1430;
    pair<double, double> BwK_0_1430;
    pair<double, double> BwK1400;
    pair<double, double> BwK_1_1270;
    pair<double, double> Bwf980_pi;
    pair<double, double> Bwf1710_k;
    pair<double, double> Bwf2200;
    pair<double, double> BwSigma;
    pair<double, double> Bwf1370_k;
    pair<double, double> BwK_0_2400;

    double Flatf980;
    double Flatf980gPiPi;
    double Flatf980gKK;

    double phaseSpace;
  };
  
  struct Psi2STo2K2PiGamEvtData
  {
    Vector4<float> psi2S_4V;
    Vector4<float>  chic0_HeliPsi2S_4V;
    Vector4<float> KpKm_HeliChic0_4V;
    Vector4<float> KpPi0_HeliChic0_4V;
    Vector4<float> KpPi1_HeliChic0_4V;
    Vector4<float> KmPi0_HeliChic0_4V;
    Vector4<float> KmPi1_HeliChic0_4V;
    Vector4<float> PiPi_HeliChic0_4V;

    Vector4<float> KpPiPi_HeliChic0_4V;
    Vector4<float> KmPiPi_HeliChic0_4V;

    Vector4<float> KKPi0_HeliChic0_4V;
    Vector4<float> KKPi1_HeliChic0_4V;

    Vector4<float> KpPi0_HeliKpPi0Pi0_4V;    
    Vector4<float> KpPi1_HeliKpPi0Pi0_4V;
    Vector4<float> KmPi0_HeliKmPi0Pi0_4V;    
    Vector4<float> KmPi1_HeliKmPi0Pi0_4V;

    Vector4<float> Kp_HeliKpPi0_4V;
    Vector4<float> Kp_HeliKpPi1_4V;
    Vector4<float> Km_HeliKmPi0_4V;
    Vector4<float> Km_HeliKmPi1_4V;

    Vector4<float> Kp_HeliKpPi0_ViaKpPiPi_4V;
    Vector4<float> Kp_HeliKpPi1_ViaKpPiPi_4V;
    Vector4<float> Km_HeliKmPi0_ViaKmPiPi_4V;
    Vector4<float> Km_HeliKmPi1_ViaKmPiPi_4V;

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfPsi; //Wigner D functions for Psi(2S) decay

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfK1400pToKstpPi1; //Wigner D functions for K+(1400)->K*+(892) pi1; K*+(892)->K+ pi0
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfK1400pToKstpPi0; //Wigner D functions for K+(1400)->K*+(892) pi0; K*+(892)->K+ pi1
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfK1400mToKstmPi1; //Wigner D functions for K-(1400)->K*-(892) pi1; K*-(892)->K- pi0
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfK1400mToKstmPi0; //Wigner D functions for K-(1400)->K*-(892) pi0; K*-(892)->K- pi1

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst1pToKpPi0; //Wigner D functions for K*1+->K+ pi0
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst1pToKpPi1; //Wigner D functions for K*1+->K+ pi1
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst1mToKmPi0; //Wigner D functions for K*1-->K- pi0
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst1mToKmPi1; //Wigner D functions for K*1-->K- pi1

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst1pToKpPi0ViaKpPiPi; //Wigner D functions for K*1+->K+ pi0
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst1pToKpPi1ViaKpPiPi; //Wigner D functions for K*1+->K+ pi1
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst1mToKmPi0ViaKmPiPi; //Wigner D functions for K*1-->K- pi0
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst1mToKmPi1ViaKmPiPi; //Wigner D functions for K*1-->K- pi1

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst2pToKpPi0; //Wigner D functions for K*2+->K+ pi0
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst2pToKpPi1; //Wigner D functions for K*2+->K+ pi1
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst2mToKmPi0; //Wigner D functions for K*2-->K- pi0
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst2mToKmPi1; //Wigner D functions for K*2-->K- pi1

//     map<Spin,map<Spin,map<Spin,complex<double> > > > DfChiToKmPi; //Wigner D functions for Chi_c1 -> (K- pi) K+
//     map<Spin,map<Spin,map<Spin,complex<double> > > > DfChiToa0Pi; //Wigner D functions for Chi_c1 ->  a(980) pi
//     map<Spin,map<Spin,map<Spin,complex<double> > > > DfChiToa2Pi; //Wigner D functions for Chi_c1 ->  a2 pi

//     map<Spin,map<Spin,map<Spin,complex<double> > > > DfKp1; //Wigner D functions for K*1+ decay
//     map<Spin,map<Spin,map<Spin,complex<double> > > > DfKp2; //Wigner D functions for K*2+ decay
//     map<Spin,map<Spin,map<Spin,complex<double> > > > DfKm1; //Wigner D functions for K*1- decay
//     map<Spin,map<Spin,map<Spin,complex<double> > > > DfKm2; //Wigner D functions for K*2- decay
//     map<Spin,map<Spin,map<Spin,complex<double> > > > Dfa2; //Wigner D functions for a2 decay
  };
  
}
#endif /* _Psi2STo2K2PiGamData_H */
