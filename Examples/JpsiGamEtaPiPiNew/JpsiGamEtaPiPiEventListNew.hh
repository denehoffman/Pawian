#ifndef _JpsiGamEtaPiPiEventListNew_H
#define _JpsiGamEtaPiPiEventListNew_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"

struct enumJpsiGamEtaPiPi4V{
  /*  enum {EtaPipPim_HeliPsi=0, PipPim_HeliPsi, EtaPip_HeliPsi, EtaPim_HeliPsi,
        gamma_HeliPsi, Eta_HeliPsi, Pip_HeliPsi, Pim_HeliPsi,
        Eta_HeliEtaPipPim, Pip_HeliEtaPipPim, Pim_HeliEtaPipPim, EtaPip_HeliEtaPipPim, EtaPim_HeliEtaPipPim, PipPim_HeliEtaPipPim,
        Pip_HeliPipPim, Pim_HeliPipPim, Pip_HeliEtaPip, Pim_HeliEtaPim, Eta_HeliEtaPip, Eta_HeliEtaPim,
        n4Vecs};*/
  enum {EtaPipPim_HeliPsi=0, PipPim_HeliPsi, EtaPip_HeliPsi, EtaPim_HeliPsi,
        gamma_HeliPsi, EtaPip_HeliEtaPipPim, EtaPim_HeliEtaPipPim, PipPim_HeliEtaPipPim,
        Pip_HeliPipPim, Pim_HeliPipPim, Pip_HeliEtaPip, Pim_HeliEtaPim, Eta_HeliEtaPip, Eta_HeliEtaPim,
        n4Vecs};

  static const std::string& name(unsigned int t)
  {
    static std::string fitName[enumJpsiGamEtaPiPi4V::n4Vecs]
      /*      ={"EtaPipPim_HeliPsi", "PipPim_HeliPsi","EtaPip_HeliPsi","EtaPim_HeliPsi",
        "gamma_HeliPsi", "Eta_HeliPsi", "Pip_HeliPsi", "Pim_HeliPsi",
        "Eta_HeliEtaPipPim", "Pip_HeliEtaPipPim", "Pim_HeliEtaPipPim", "EtaPip_HeliEtaPipPim", "EtaPim_HeliEtaPipPim", "PipPim_HeliEtaPipPim",
        "Pip_HeliPipPim", "Pim_HeliPipPim", "Pip_HeliEtaPip", "Pim_HeliEtaPim", "Eta_HeliEtaPip", "Eta_HeliEtaPim",
	};*/
    ={"EtaPipPim_HeliPsi", "PipPim_HeliPsi","EtaPip_HeliPsi","EtaPim_HeliPsi",
      "gamma_HeliPsi", "EtaPip_HeliEtaPipPim", "EtaPim_HeliEtaPipPim","PipPim_HeliEtaPipPim",
      "Pip_HeliPipPim", "Pim_HeliPipPim", "Pip_HeliEtaPip", "Pim_HeliEtaPim", "Eta_HeliEtaPip", "Eta_HeliEtaPim"
    };

    
    if (t<0 || t>=enumJpsiGamEtaPiPi4V::n4Vecs) assert(0);
    return fitName[t];
  }
};

struct enumJpsiGamEtaPiPiDfunc{
  enum {etapipi=0, XToAplusPiminus, XToAminusPiplus, XTofEta,
        Aplus, Aminus, fToPiPi, nDfuncts};

  static const std::string& name(unsigned int t)
  {
    static std::string fitName[enumJpsiGamEtaPiPiDfunc::nDfuncts]
      ={"etapipi", "XToAplusPiminus", "XToAminusPiplus", "XTofEta", 
        "Aplus", "Aminus", "fToPiPi"};
    if (t<0 || t>=enumJpsiGamEtaPiPiDfunc::nDfuncts) assert(0);
    return fitName[t];
  }
};



class EventList;

class JpsiGamEtaPiPiEventListNew : public EvtDataBaseListNew {

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamEtaPiPiEventListNew();



  /** Destructor */
  virtual ~JpsiGamEtaPiPiEventListNew();

  // Getters:
  

protected:


  virtual void read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList, double& evtWeightSum, int maxEvts);
  
private:


};

#endif 
