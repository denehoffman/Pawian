#ifndef _JpsiToOmegaPhiGamEventList_H
#define _JpsiToOmegaPhiGamEventList_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"

struct enumJpsiGamX4V{
  enum { Psi=0,  V4_PipPimPi0KpKm_HeliPsi, V4_PipPimPi0_HeliPsi, V4_KpKm_HeliPsi, V4_gamma_HeliPsi,
	 V4_Pip_HeliPsi, V4_Pim_HeliPsi, V4_Pi0_HeliPsi, V4_Kp_HeliPsi, V4_Km_HeliPsi,
	 V4_PipPimPi0_HeliPipPimPi0KpKm, V4_KpKm_HeliPipPimPi0KpKm, V4_Pip_HeliPipPimPi0, V4_Kp_HeliKpKm, n4Vecs};
  
  static const std::string& name(unsigned int t)
  {
    static std::string fitName[enumJpsiGamX4V::n4Vecs]
      ={ "Psi",  "PipPimPi0KpKm_HeliPsi", "PipPimPi0_HeliPsi", "KpKm_HeliPsi", "gamma_HeliPsi",
	 "Pip_HeliPsi", "Pim_HeliPsi", "Pi0_HeliPsi", "Kp_HeliPsi", "Km_HeliPsi",
	 "PipPimPi0_HeliPipPimPi0KpKm", "KpKm_HeliPipPimPi0KpKm", "Pip_HeliPipPimPi0", "Kp_HeliKpKm"  };
    
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

class JpsiToOmegaPhiGamEventList : public EvtDataBaseListNew {

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiToOmegaPhiGamEventList();



  /** Destructor */
  virtual ~JpsiToOmegaPhiGamEventList();

  // Getters:
  

protected:


  virtual void read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList, double& evtWeightSum, int maxEvts);
  
private:


};

#endif 
