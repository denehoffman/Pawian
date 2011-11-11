#ifndef _Psi2SToKpKmPiGamData_H
#define _Psi2SToKpKmPiGamData_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"

struct enumChic1ToKpKmPiGamData{
  enum {Psi2S=0, Chic1_HeliPsi2S, KpKm_HeliChic1, KpPi_HeliChic1, KmPi_HeliChic1, KpKmPi0_HeliChic1,
	Kp_HeliKpKm, Kp_HeliKpPi, Km_HeliKmPi, n4Vecs,
	dfPsi=n4Vecs, DfChiToKpPi, DfChiToKmPi, DfChiToa0Pi, DfChiToa2Pi, DfChiToKKPi, 
	DfKp1, DfKp2, DfKm1, DfKm2, Dfa2, nDfuncts};

  static const std::string& name(unsigned int t)
  {
    static std::string fitName[enumChic1ToKpKmPiGamData::nDfuncts]
      ={"Psi2S", "Chic1_HeliPsi2S", "KpKm_HeliChic1", "KpPi_HeliChic1", "KmPi_HeliChic1", "KpKmPi0_HeliChic1",
	"Kp_HeliKpKm", "Kp_HeliKpPi", "Km_HeliKmPi",
	"dfPsi", "DfChiToKpPi", "DfChiToKmPi", "DfChiToa0Pi", "DfChiToa2Pi", "DfChiToKKPi", 
	"DfKp1", "DfKp2", "DfKm1", "DfKm2", "Dfa2"};
    if (t<0 || t>=enumChic1ToKpKmPiGamData::nDfuncts) assert(0);
    return fitName[t];
  }
};



namespace Psi2SToKpKmPiGamData {

  struct Psi2SToKpKmPiGamEvtData
  {
    Vector4<float> psi2S_4V;
    Vector4<float>  chic1_HeliPsi2S_4V;
    Vector4<float> KpKm_HeliChic1_4V;
    Vector4<float> KpPi_HeliChic1_4V;
    Vector4<float> KmPi_HeliChic1_4V;
    Vector4<float> KpKmPi0_HeliChic1_4V;
    Vector4<float> Kp_HeliKpKm_4V;
    Vector4<float> Kp_HeliKpPi_4V;
    Vector4<float> Km_HeliKmPi_4V;
    map<Spin,map<Spin,map<Spin,complex<double> > > > dfPsi; //Wigner D functions for Psi(2S) decay
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfChiToKpPi; //Wigner D functions for Chi_c1 -> (K+ pi) K-
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfChiToKmPi; //Wigner D functions for Chi_c1 -> (K- pi) K+
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfChiToa0Pi; //Wigner D functions for Chi_c1 ->  a(980) pi
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfChiToa2Pi; //Wigner D functions for Chi_c1 ->  a2 pi
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfChiToKKPi; //Wigner D functions for Chi_c1 ->  K+ K- pi

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKp1; //Wigner D functions for K*1+ decay
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKp2; //Wigner D functions for K*2+ decay
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKm1; //Wigner D functions for K*1- decay
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKm2; //Wigner D functions for K*2- decay
    map<Spin,map<Spin,map<Spin,complex<double> > > > Dfa2; //Wigner D functions for a2 decay
  };
  
}
#endif /* _Psi2SToKpKmPiGamData_H */
