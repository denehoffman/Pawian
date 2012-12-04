#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/D0ToKsPipPim/D0ToKsPipPimEventList.hh"

#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


D0ToKsPipPimEventList::D0ToKsPipPimEventList() :
  EvtDataBaseListNew()
{
}


D0ToKsPipPimEventList::~D0ToKsPipPimEventList()
{
}

void D0ToKsPipPimEventList::read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList, double& evtWeightSum, int maxEvts){
  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = evtList.nextEvent())){
    if (evtCount>= maxEvts) break;
    if (evtCount%10000 == 0) Info << "4vec calculation for event " << evtCount ;  // << endmsg;

    Vector4<float> ks_lab = *(anEvent->p4(0));
    Vector4<float> pip_lab = *(anEvent->p4(1));
    Vector4<float> pim_lab = *(anEvent->p4(2));
    
    Vector4<float> V4_D0_lab = ks_lab+pip_lab+pim_lab;
    Vector4<float> V4_pippim_lab = pip_lab+pim_lab;
    Vector4<float> V4_KsPip_lab = ks_lab+pip_lab;
    Vector4<float> V4_KsPim_lab = ks_lab+pim_lab;

    if (evtCount%10000 == 0){
      Info << "D0 4vec" << "\n"
	   << " px: " << V4_D0_lab.Px() <<"\t"
	   << " py: " << V4_D0_lab.Py() <<"\t"
	   << " pz: " << V4_D0_lab.Pz() <<"\t"
	   << " e : " << V4_D0_lab.E() << "\t"
	   << " m : " << V4_D0_lab.M() ;  // << endmsg;
    }

    Vector4<float> V4_D0_Heli(V4_D0_lab);
    V4_D0_Heli.Boost(V4_D0_lab);
    
    Vector4<float>  V4_Ks_HeliD0( ks_lab   );
    V4_Ks_HeliD0.Boost(V4_D0_lab);

    Vector4<float>  V4_pip_HeliD0( pip_lab   );
    V4_pip_HeliD0.Boost(V4_D0_lab); 

    Vector4<float>  V4_pim_HeliD0( pip_lab   );
    V4_pim_HeliD0.Boost(V4_D0_lab);          

    Vector4<float>  V4_KsPip_HeliD0( ks_lab+pip_lab   );
    V4_KsPip_HeliD0.Boost(V4_D0_lab);

    Vector4<float>  V4_KsPim_HeliD0( ks_lab+pim_lab   );
    V4_KsPim_HeliD0.Boost(V4_D0_lab);

    Vector4<float>  V4_PipPim_HeliD0( pip_lab+pim_lab   );
    V4_PipPim_HeliD0.Boost(V4_D0_lab);

    Vector4<float> V4_Pip_HeliPipPim=helicityVec(V4_D0_lab, V4_pippim_lab, pip_lab);
    Vector4<float> V4_Pip_HeliKsPip=helicityVec(V4_D0_lab, V4_KsPip_lab, pip_lab);
    Vector4<float> V4_Pim_HeliKsPim=helicityVec(V4_D0_lab, V4_KsPim_lab, pim_lab);
  
    EvtDataNew* evtData=new EvtDataNew();
    evtData->FourVecsDec[enumD0KsPiPi4V::V4_D0_Heli] = V4_D0_Heli;
    evtData->FourVecsDec[enumD0KsPiPi4V::V4_Ks_HeliD0] = V4_Ks_HeliD0;
    evtData->FourVecsDec[enumD0KsPiPi4V::V4_Pip_HeliD0] = V4_pip_HeliD0;
    evtData->FourVecsDec[enumD0KsPiPi4V::V4_Pim_HeliD0] = V4_pim_HeliD0;
    evtData->FourVecsDec[enumD0KsPiPi4V::V4_PipPim_HeliD0] = V4_PipPim_HeliD0;
    evtData->FourVecsDec[enumD0KsPiPi4V::V4_KsPip_HeliD0] = V4_KsPip_HeliD0;
    evtData->FourVecsDec[enumD0KsPiPi4V::V4_KsPim_HeliD0] = V4_KsPim_HeliD0;
    evtData->FourVecsDec[enumD0KsPiPi4V::V4_Pip_HeliPipPim] = V4_Pip_HeliPipPim;
    evtData->FourVecsDec[enumD0KsPiPi4V::V4_Pip_HeliKsPip] = V4_Pip_HeliKsPip;
    evtData->FourVecsDec[enumD0KsPiPi4V::V4_Pim_HeliKsPim] = V4_Pim_HeliKsPim;
 
    // calculate and store Wignerd functions for X -> pi+ pi-, X -> Ks pi+ and X -> Ks pi- 

    for (Spin jX=0; jX<2; jX++){
      evtData->WignerDsDec[enumD0KsPiPiDfunc::df_XToPipPim][jX][0][0]=Wigner_D(0.,V4_Pip_HeliPipPim.Theta(),0,jX,0,0);
      evtData->WignerDsDec[enumD0KsPiPiDfunc::df_XToKsPip][jX][0][0]=Wigner_D(0.,V4_Pip_HeliKsPip.Theta(),0,jX,0,0);
      evtData->WignerDsDec[enumD0KsPiPiDfunc::df_XToKsPip][jX][0][0]=Wigner_D(0.,V4_Pim_HeliKsPim.Theta(),0,jX,0,0);
    }

    evtData->evtWeight=anEvent->Weight();
    evtData->evtNo=_evtNoAll;
    theEvtList.push_back(evtData);
    
    evtWeightSum += anEvent->Weight();    
    ++evtCount;
    ++_evtNoAll;
  }
}


