#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiToPhiPhiGam/JpsiToPhiPhiGamEventList.hh"

#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


JpsiToPhiPhiGamEventList::JpsiToPhiPhiGamEventList() :
  EvtDataBaseListNew()
{
}


JpsiToPhiPhiGamEventList::~JpsiToPhiPhiGamEventList()
{
}

void JpsiToPhiPhiGamEventList::read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList, double& evtWeightSum, int maxEvts){
  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = evtList.nextEvent())){
    if (evtCount>= maxEvts) break;
    if (evtCount%10000 == 0) Info << "4vec calculation for event " << evtCount ;  // << endmsg;

    Vector4<float> gam = *(anEvent->p4(0));
    Vector4<float> ks = *(anEvent->p4(1));
    Vector4<float> kl = *(anEvent->p4(2));
    Vector4<float> kp = *(anEvent->p4(3));
    Vector4<float> km = *(anEvent->p4(4));
    
    Vector4<float> V4_psi = gam+ks+kl+kp+km;

    if (evtCount%10000 == 0){
      Info << "psi 4vec" << "\n"
	   << " px: " << V4_psi.Px() <<"\t"
	   << " py: " << V4_psi.Py() <<"\t"
	   << " pz: " << V4_psi.Pz() <<"\t"
	   << " e : " << V4_psi.E() << "\t"
	   << " m : " << V4_psi.M() ;  // << endmsg;
    }
    
     Vector4<float>  V4_all_Lab( ks+kl+kp+km+gam   );
     Vector4<float>  V4_KsKlKpKm_Lab( ks+kl+kp+km   );
     Vector4<float>  V4_KsKl_Lab( ks+kl   );
     Vector4<float>  V4_Ks_Lab( ks   );
     Vector4<float>  V4_KpKm_Lab( kp+km   );
     Vector4<float>  V4_Kp_Lab( kp   );
     Vector4<float>  V4_Km_Lab( km   );
     Vector4<float>  V4_Kl_Lab( kl   );
          
     Vector4<float>  V4_KsKlKpKm_HeliPsi( ks+kl+kp+km   );
     V4_KsKlKpKm_HeliPsi.Boost(V4_psi);
     
     Vector4<float>  V4_gamma_HeliPsi( gam   );
     V4_gamma_HeliPsi.Boost(V4_psi);
     Vector4<float>  V4_KpKm_HeliPsi( kp+km  );
     V4_KpKm_HeliPsi.Boost(V4_psi);
     Vector4<float>  V4_KsKl_HeliPsi( ks+kl  );
     V4_KsKl_HeliPsi.Boost(V4_psi);

     Vector4<float>  V4_Kl_HeliPsi( kl   );
     V4_Kl_HeliPsi.Boost( V4_psi );
     Vector4<float>  V4_Ks_HeliPsi( ks   );
     V4_Ks_HeliPsi.Boost( V4_psi );
     Vector4<float>  V4_Kp_HeliPsi( kp   );
     V4_Kp_HeliPsi.Boost( V4_psi );
     Vector4<float>  V4_Km_HeliPsi( km   );
     V4_Km_HeliPsi.Boost( V4_psi );

     Vector4<float>  V4_KsKl_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_KsKl_Lab);
     Vector4<float>  V4_KpKm_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_KpKm_Lab);     
     Vector4<float>  V4_Ks_HeliKsKl=helicityVec(V4_KsKlKpKm_Lab, V4_KsKl_Lab, V4_Ks_Lab);     
     Vector4<float>  V4_Kp_HeliKpKm=helicityVec(V4_KsKlKpKm_Lab, V4_KpKm_Lab, V4_Kp_Lab);

     Vector4<float>  V4_Kp_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_Kp_Lab);
     Vector4<float>  V4_Km_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_Km_Lab);
     
     Vector4<float> V4_normKpKmDecHeliKsKlKpKm
       (0.5*(V4_Kp_HeliKsKlKpKm.T()+V4_Km_HeliKsKlKpKm.T()),
	V4_Km_HeliKsKlKpKm.Y()*V4_Kp_HeliKsKlKpKm.Z()-V4_Km_HeliKsKlKpKm.Z()*V4_Kp_HeliKsKlKpKm.Y(),
	V4_Km_HeliKsKlKpKm.Z()*V4_Kp_HeliKsKlKpKm.X()-V4_Km_HeliKsKlKpKm.X()*V4_Kp_HeliKsKlKpKm.Z(),
	V4_Km_HeliKsKlKpKm.X()*V4_Kp_HeliKsKlKpKm.Y()-V4_Km_HeliKsKlKpKm.Y()*V4_Kp_HeliKsKlKpKm.X());
     
     Vector4<float>  V4_Kl_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_Kl_Lab);
     Vector4<float>  V4_Ks_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_Ks_Lab);
     
     Vector4<float> V4_normKsKlDecHeliKsKlKpKm
       (0.5*(V4_Kl_HeliKsKlKpKm.T()+V4_Ks_HeliKsKlKpKm.T()),
	V4_Ks_HeliKsKlKpKm.Y()*V4_Kl_HeliKsKlKpKm.Z()-V4_Ks_HeliKsKlKpKm.Z()*V4_Kl_HeliKsKlKpKm.Y(),
	V4_Ks_HeliKsKlKpKm.Z()*V4_Kl_HeliKsKlKpKm.X()-V4_Ks_HeliKsKlKpKm.X()*V4_Kl_HeliKsKlKpKm.Z(),
	V4_Ks_HeliKsKlKpKm.X()*V4_Kl_HeliKsKlKpKm.Y()-V4_Ks_HeliKsKlKpKm.Y()*V4_Kl_HeliKsKlKpKm.X());
     
 
     EvtDataNew* evtData=new EvtDataNew();
     evtData->FourVecsProd[enumProd4V::Psi] = V4_psi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_KsKlKpKm_HeliPsi]=V4_KsKlKpKm_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_KsKl_HeliPsi]=V4_KsKl_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_KpKm_HeliPsi]=V4_KpKm_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_gamma_HeliPsi]=V4_gamma_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Ks_HeliPsi]=V4_Ks_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Kl_HeliPsi]=V4_Kl_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Kp_HeliPsi]=V4_Kp_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Km_HeliPsi]=V4_Km_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_KsKl_HeliKsKlKpKm]=V4_KsKl_HeliKsKlKpKm;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_KpKm_HeliKsKlKpKm]=V4_KpKm_HeliKsKlKpKm;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Ks_HeliKsKl]=V4_Ks_HeliKsKl;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Kp_HeliKpKm]= V4_Kp_HeliKpKm;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_normKpKmDecHeliKsKlKpKm]=V4_normKpKmDecHeliKsKlKpKm;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_normKsKlDecHeliKsKlKpKm]=V4_normKsKlDecHeliKsKlKpKm;

    
    // calculate and store WignerD functions for Psi -> X gamma (JPC =0-+, 0++, 2++)

    Spin jPsi=1;
    for (Spin M=-1; M<=1; M=M+2){
      for (Spin lam=-1; lam<=1; lam++){
	evtData->WignerDsProd[enumProdDfunc::Psi][jPsi][M][lam]=Wigner_D(0.,V4_KsKlKpKm_HeliPsi.Theta(),0,jPsi,M,lam);
      }
    }

  //WignerD functions for X -> phi phi

    for (Spin J_X=0; J_X<=2; J_X++){
      for (Spin lam_X=-J_X; lam_X<=J_X; lam_X++){
	for (Spin lamPhi1mlamPhi2=-J_X; lamPhi1mlamPhi2<=J_X; lamPhi1mlamPhi2++){
	evtData->WignerDsDec[enumJpsiGamXDfunc::Df_XToPhiPhi][J_X][lam_X][lamPhi1mlamPhi2]
	  =Wigner_D(V4_KsKl_HeliKsKlKpKm.Phi(),V4_KsKl_HeliKsKlKpKm.Theta(),0,J_X,lam_X,lamPhi1mlamPhi2);
      }
    }
    }
    
    //WignerD function for phi -> K+ K- and phi -> KS KL
    Spin phiSpin=1;
    for(Spin M=-phiSpin; M<=phiSpin; M++){
      Spin lam=0;
      evtData->WignerDsDec[enumJpsiGamXDfunc::Df_PhiToKsKl][phiSpin][M][lam] = Wigner_D(V4_Ks_HeliKsKl.Phi(),V4_Ks_HeliKsKl.Theta(), 0,phiSpin,M,lam);
      evtData->WignerDsDec[enumJpsiGamXDfunc::Df_PhiToKpKm][phiSpin][M][lam] = Wigner_D(V4_Kp_HeliKpKm.Phi(),V4_Kp_HeliKpKm.Theta(), 0,phiSpin,M,lam);
     }    
    
    evtData->evtWeight=anEvent->Weight();
    theEvtList.push_back(evtData);
    
    evtWeightSum += anEvent->Weight();    
    ++evtCount;
  }
}


