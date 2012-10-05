#ifndef _JpsiToPhiPhiGamEventList_H
#define _JpsiToPhiPhiGamEventList_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"

struct enumJpsiGamX4V{
  enum {V4_KsKlKpKm_HeliPsi=0, V4_KsKl_HeliPsi,V4_KpKm_HeliPsi, V4_gamma_HeliPsi,
        V4_Ks_HeliPsi, V4_Kl_HeliPsi, V4_Kp_HeliPsi, V4_Km_HeliPsi,
        V4_KsKl_HeliKsKlKpKm,V4_KpKm_HeliKsKlKpKm ,V4_Ks_HeliKsKl,V4_Kp_HeliKpKm,
        V4_normKpKmDecHeliKsKlKpKm, V4_normKsKlDecHeliKsKlKpKm, n4Vecs};

  static const std::string& name(unsigned int t)
  {
    static std::string fitName[enumJpsiGamX4V::n4Vecs]
      ={"KsKlKpKm_HeliPsi", "KsKl_HeliPsi,KpKm_HeliPsi", "gamma_HeliPsi",
        "Ks_HeliPsi, Kl_HeliPsi", "Kp_HeliPsi", "Km_HeliPsi",
        "KsKl_HeliKsKlKpKm","KpKm_HeliKsKlKpKm" ,"Ks_HeliKsKl","Kp_HeliKpKm",
        "normKpKmDecHeliKsKlKpKm", "normKsKlDecHeliKsKlKpKm"
    };

    if (t<0 || t>=enumJpsiGamX4V::n4Vecs) assert(0);
    return fitName[t];
  }
};

struct enumJpsiGamXDfunc{
  enum {Df_XToPhiPhi1=0, Df_XToPhiPhi2, Df_PhiToKsKl,Df_PhiToKpKm, nDfuncts};

  static const std::string& name(unsigned int t)
  {
    static std::string fitName[enumJpsiGamXDfunc::nDfuncts]
      ={"XToPhiPhi1", "XToPhiPhi2", "PhiToKsKl", "PhiToKpKm"};
    if (t<0 || t>=enumJpsiGamXDfunc::nDfuncts) assert(0);
    return fitName[t];
  }
};



class EventList;

class JpsiToPhiPhiGamEventList : public EvtDataBaseListNew {

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiToPhiPhiGamEventList();



  /** Destructor */
  virtual ~JpsiToPhiPhiGamEventList();

  // Getters:
  

protected:


  virtual void read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList, double& evtWeightSum, int maxEvts);
  
private:


};

#endif 
