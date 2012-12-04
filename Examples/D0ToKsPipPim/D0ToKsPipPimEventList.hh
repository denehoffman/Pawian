#pragma once

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"

struct enumD0KsPiPi4V{
  enum {V4_D0_Heli=0,  V4_Ks_HeliD0, V4_Pip_HeliD0, V4_Pim_HeliD0,
	V4_PipPim_HeliD0, V4_KsPip_HeliD0, V4_KsPim_HeliD0,
	V4_Pip_HeliPipPim, V4_Pip_HeliKsPip, V4_Pim_HeliKsPim, 
        n4Vecs};

  static const std::string& name(unsigned int t)
  {
    static std::string fitName[enumD0KsPiPi4V::n4Vecs]
      ={"D0_HeliD0", "Ks_HeliD0", "Pip_HeliD0", "Pim_HeliD0",
	"PipPim_HeliD0", "KsPip_HeliD0", "KsPim_HeliD0",
	"Pip_HeliPipPim", "Pip_HeliKsPip","Pim_HeliKsPim"};

    if (t<0 || t>=enumD0KsPiPi4V::n4Vecs) assert(0);
    return fitName[t];
  }
};

struct enumD0KsPiPiDfunc{
  enum {df_XToPipPim=0, df_XToKsPip, df_XToKsPim, nDfuncts};

  static const std::string& name(unsigned int t)
  {
    static std::string fitName[enumD0KsPiPiDfunc::nDfuncts]
      ={"XToPipPim", "XToKsPip", "XToKsPim"};
    if (t<0 || t>=enumD0KsPiPiDfunc::nDfuncts) assert(0);
    return fitName[t];
  }
};



class EventList;

class D0ToKsPipPimEventList : public EvtDataBaseListNew {

public:

  // create/copy/destroy:

  ///Constructor 
  D0ToKsPipPimEventList();



  /** Destructor */
  virtual ~D0ToKsPipPimEventList();

  // Getters:
  

protected:


  virtual void read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList, double& evtWeightSum, int maxEvts);
  
private:


};

