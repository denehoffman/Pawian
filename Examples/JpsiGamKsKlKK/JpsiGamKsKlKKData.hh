#ifndef _JpsiGamKsKlKKData_H
#define _JpsiGamKsKlKKData_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"

namespace JpsiGamKsKlKKData {

  struct JpsiGamKsKlKKEvtData
  {
    

    Vector4<float> V4_Psi;
    
    Vector4<float> V4_KsKlKpKm_HeliPsi;
    Vector4<float> V4_KsKl_HeliPsi;
    Vector4<float> V4_KpKm_HeliPsi;
    Vector4<float> V4_gamma_HeliPsi;
    
    Vector4<float> V4_KsKl_HeliKsKlKpKm;
    Vector4<float> V4_KpKm_HeliKsKlKpKm;
    Vector4<float> V4_Ks_HeliKsKl;
    Vector4<float> V4_Kp_HeliKpKm;


    map<Spin,map<Spin,map<Spin,complex<double> > > > df_Psi;

    map<Spin,map<Spin,map<Spin,complex<double> > > > Df_Spin0;
    map<Spin,map<Spin,map<Spin,complex<double> > > > Df_Spin2;
    map<Spin,map<Spin,map<Spin,complex<double> > > > Df_KsKl;
    map<Spin,map<Spin,map<Spin,complex<double> > > > Df_KpKm;
    
    

  };
  
}
#endif /* _JpsiGamKsKlKKData_H */
