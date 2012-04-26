#ifndef _JpsiGamEtaPiPiData_H
#define _JpsiGamEtaPiPiData_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"


struct enumJpsiGamEtaPiPiData{
  enum {V4_Psi=0, V4_EtaPipPim_HeliPsi, V4_PipPim_HeliPsi, V4_EtaPip_HeliPsi, V4_EtaPim_HeliPsi,
	V4_gamma_HeliPsi, V4_Eta_HeliPsi, V4_Pip_HeliPsi, V4_Pim_HeliPsi,
	V4_Eta_HeliEtaPipPim, V4_Pip_HeliEtaPipPim, V4_Pim_HeliEtaPipPim, V4_EtaPip_HeliEtaPipPim, V4_EtaPim_HeliEtaPipPim, V4_PipPim_HeliEtaPipPim,
	V4_Pip_HeliPipPim, V4_Pim_HeliPipPim, V4_Pip_HeliEtaPip, V4_Pim_HeliEtaPim, V4_Eta_HeliEtaPip, V4_Eta_HeliEtaPim,
	V4_normPipPimDecHeliEtaPipPim, V4_normEtaPipDecHeliEtaPipPim,
	n4Vecs, Df_Psi=n4Vecs, Df_etapipidec, Df_AplusPiminusdec, Df_AminusPiplusdec, nDfuncts
  };
  
  static const std::string& name(unsigned int t)
  {
    static std::string fitName[enumJpsiGamEtaPiPiData::nDfuncts]
      ={"Psi", "EtaPipPim_HeliPsi", "PipPim_HeliPsi","EtaPip_HeliPsi","EtaPim_HeliPsi",
    	"gamma_HeliPsi", "Eta_HeliPsi", "Pip_HeliPsi", "Pim_HeliPsi",
    	"Eta_HeliEtaPipPim", "Pip_HeliEtaPipPim", "Pim_HeliEtaPipPim", "EtaPip_HeliEtaPipPim", "EtaPim_HeliEtaPipPim", "PipPim_HeliEtaPipPim",
    	"Pip_HeliPipPim", "Pim_HeliPipPim", "Pip_HeliEtaPip", "Pim_HeliEtaPim", "Eta_HeliEtaPip", "Eta_HeliEtaPim",
    	"Df_Psi","Df_etapipidec", "Df_apidec", "Df_AplusPiminusdec", "Df_AminusPiplusdec"
    };
    if (t<0 || t>=enumJpsiGamEtaPiPiData::nDfuncts) assert(0);
    return fitName[t];
  }
};

#endif /* _JpsiGamEtaPiPiData_H */
