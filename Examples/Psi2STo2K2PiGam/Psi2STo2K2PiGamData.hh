#ifndef _Psi2STo2K2PiGamData_H
#define _Psi2STo2K2PiGamData_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"


namespace Psi2STo2K2PiGamData {

  struct Psi2STo2K2PiGamEvtData
  {
    int evtNo;
    Vector4<float> psi2S_4V;
    Vector4<float> chic0_HeliPsi2S_4V;
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

    Vector4<float> PiPi_HeliKpPi0Pi0_4V;    
    Vector4<float> PiPi_HeliKmPi0Pi0_4V;

    Vector4<float> Pi0_HeliPi0Pi0_ViaKpPi0Pi0_4V;
    Vector4<float> Pi0_HeliPi0Pi0_ViaKmPi0Pi0_4V;

    Vector4<float> Kp_HeliKpPi0_4V;
    Vector4<float> Kp_HeliKpPi1_4V;
    Vector4<float> Km_HeliKmPi0_4V;
    Vector4<float> Km_HeliKmPi1_4V;

    Vector4<float> Km_HeliKmKp_4V;
    Vector4<float> Pi0_HeliPi0Pi0_4V;

    Vector4<float> Kp_HeliKpPi0_ViaKpPiPi_4V;
    Vector4<float> Kp_HeliKpPi1_ViaKpPiPi_4V;
    Vector4<float> Km_HeliKmPi0_ViaKmPiPi_4V;
    Vector4<float> Km_HeliKmPi1_ViaKmPiPi_4V;

    Vector4<float> KK_HeliKKPi0_4V;
    Vector4<float> KK_HeliKKPi1_4V;
 
    Vector4<float> KpPi0_HeliKKPi0_4V;    
    Vector4<float> KpPi1_HeliKKPi1_4V;
    Vector4<float> KmPi0_HeliKKPi0_4V;    
    Vector4<float> KmPi1_HeliKKPi1_4V;

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

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfPi2Tof2Pi0; //Wigner D functions 
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfPi2Tof2Pi1; //Wigner D functions 

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfPi2ToKstarpK_pi0; //Wigner D functions 
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfPi2ToKstarpK_pi1; //Wigner D functions 
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfPi2ToKstarmK_pi0; //Wigner D functions 
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfPi2ToKstarmK_pi1; //Wigner D functions 

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKp2ToKppi0; //Wigner D functions 
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKp2ToKppi1; //Wigner D functions
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKm2ToKmpi0; //Wigner D functions 
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKm2ToKmpi1; //Wigner D functions

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfK2pToKpPi0ViaKpPiPi;
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfK2pToKpPi1ViaKpPiPi;
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfK2mToKmPi0ViaKmPiPi;
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfK2mToKmPi1ViaKmPiPi;

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst1pToKpPi0ViaKKPi0;
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst1pToKpPi1ViaKKPi1;
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst1mToKmPi0ViaKKPi0;
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKst1mToKmPi1ViaKKPi1;

    map<Spin,map<Spin,map<Spin,complex<double> > > > Dff2ToKK; //Wigner D functions for f_2->K+ K-
    map<Spin,map<Spin,map<Spin,complex<double> > > > Dff2ToPiPi; //Wigner D functions for f_2->pi pi

    map<Spin,map<Spin,map<Spin,complex<double> > > > Dff2ToKKviaKKpi0;
    map<Spin,map<Spin,map<Spin,complex<double> > > > Dff2ToKKviaKKpi1;

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfK1pTof0Kp;
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfK1mTof0Km;

    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKjpTof0Kp;
    map<Spin,map<Spin,map<Spin,complex<double> > > > DfKjmTof0Km;

    bool operator<(const Psi2STo2K2PiGamEvtData& compare) const{
      bool result=false;
      if ( evtNo < compare.evtNo) result=true;
      return result;
    }  


  };
  
}
#endif /* _Psi2STo2K2PiGamData_H */
