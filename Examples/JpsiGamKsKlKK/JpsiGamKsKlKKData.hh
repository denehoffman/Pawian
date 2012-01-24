#ifndef _JpsiGamKsKlKKData_H
#define _JpsiGamKsKlKKData_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"


struct enumJpsiGamKsKlKKData{
  enum {V4_Psi=0, V4_KsKlKpKm_HeliPsi, V4_KsKl_HeliPsi,V4_KpKm_HeliPsi, V4_gamma_HeliPsi,
	V4_Ks_HeliPsi, V4_Kl_HeliPsi, V4_Kp_HeliPsi, V4_Km_HeliPsi,
	V4_KsKl_HeliKsKlKpKm,V4_KpKm_HeliKsKlKpKm ,V4_Ks_HeliKsKl,V4_Kp_HeliKpKm,
	V4_normKpKmDecHeliKsKlKpKm, V4_normKsKlDecHeliKsKlKpKm,
	n4Vecs, Df_Psi=n4Vecs, Df_Spin0,Df_Spin2, Df_KsKl,Df_KpKm, nDfuncts};
  
  static const std::string& name(unsigned int t)
  {
    static std::string fitName[enumJpsiGamKsKlKKData::nDfuncts]
      ={"Psi", "KsKlKpKm_HeliPsi", "KsKl_HeliPsi","KpKm_HeliPsi","gamma_HeliPsi",
	"Ks_HeliPsi", "Kl_HeliPsi", "Kp_HeliPsi", "Km_HeliPsi",
	"KsKl_HeliKsKlKpKm","KpKm_HeliKsKlKpKm" ,"Ks_HeliKsKl","Kp_HeliKpKm",
	"normKpKmDecHeliKsKlKpKm", "normKsKlDecHeliKsKlKpKm",
	"Df_Psi","Df_Spin0","Df_Spin2", "Df_KsKl","Df_KpKm"};
    if (t<0 || t>=enumJpsiGamKsKlKKData::nDfuncts) assert(0);
    return fitName[t];
  }
};

#endif /* _JpsiGamKsKlKKData_H */
