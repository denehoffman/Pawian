#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiEventListNew.hh"

#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


JpsiGamEtaPiPiEventListNew::JpsiGamEtaPiPiEventListNew(EventList& evtListData, EventList& evtListMc)
{
  read4Vecs(evtListData, _evtDataList);
  read4Vecs(evtListMc, _mcDataList);
}


JpsiGamEtaPiPiEventListNew::~JpsiGamEtaPiPiEventListNew()
{
}

void JpsiGamEtaPiPiEventListNew::read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList){
  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = evtList.nextEvent())){
    if (evtCount%10000 == 0) Info << "4vec calculation for event " << evtCount ;  // << endmsg;
    
    Vector4<float> gam = *(anEvent->p4(0));
    Vector4<float> eta = *(anEvent->p4(1));
    Vector4<float> pip = *(anEvent->p4(2));
    Vector4<float> pim = *(anEvent->p4(3));
    
    Vector4<float> V4_psi = gam+eta+pip+pim;
    
    if (evtCount%10000 == 0){
      Info << "psi 4vec" << "\n"
	   << " px: " << V4_psi.Px() <<"\t"
	   << " py: " << V4_psi.Py() <<"\t"
	   << " pz: " << V4_psi.Pz() <<"\t"
	   << " e : " << V4_psi.E() << "\t"
	   << " m : " << V4_psi.M() ;  // << endmsg;
    }
    
    Vector4<float>  V4_all_Lab(gam+eta+pip+pim);
    Vector4<float>  V4_EtaPipPim_Lab(eta+pip+pim);
    Vector4<float>  V4_EtaPip_Lab(eta+pip);
    Vector4<float>  V4_EtaPim_Lab(eta+pim);
    Vector4<float>  V4_PipPim_Lab(pip+pim);
    Vector4<float>  V4_Pip_Lab(pip);
    Vector4<float>  V4_Pim_Lab(pim);
    Vector4<float>  V4_Eta_Lab(eta);
    
    Vector4<float>  V4_EtaPipPim_HeliPsi(eta+pip+pim);
    V4_EtaPipPim_HeliPsi.Boost(V4_psi);
    
    Vector4<float>  V4_gamma_HeliPsi(gam);
    V4_gamma_HeliPsi.Boost(V4_psi);
    Vector4<float>  V4_PipPim_HeliPsi(pip+pim);
    V4_PipPim_HeliPsi.Boost(V4_psi);
    Vector4<float>  V4_EtaPip_HeliPsi(eta+pip);
    V4_EtaPip_HeliPsi.Boost(V4_psi);
    Vector4<float>  V4_EtaPim_HeliPsi(eta+pim);
    V4_EtaPim_HeliPsi.Boost(V4_psi);
    
    Vector4<float>  V4_Eta_HeliPsi(eta);
    V4_Eta_HeliPsi.Boost(V4_psi);
    Vector4<float>  V4_Pip_HeliPsi(pip);
    V4_Pip_HeliPsi.Boost(V4_psi);
    Vector4<float>  V4_Pim_HeliPsi(pim);
    V4_Pim_HeliPsi.Boost(V4_psi);
    
    
    Vector4<float>  V4_Eta_HeliEtaPipPim = helicityVec(V4_all_Lab, V4_EtaPipPim_Lab, V4_Eta_Lab);
    Vector4<float>  V4_Pip_HeliEtaPipPim = helicityVec(V4_all_Lab, V4_EtaPipPim_Lab, V4_Pip_Lab);
    Vector4<float>  V4_Pim_HeliEtaPipPim = helicityVec(V4_all_Lab, V4_EtaPipPim_Lab, V4_Pim_Lab);
    
    Vector4<float>  V4_EtaPip_HeliEtaPipPim = helicityVec(V4_all_Lab, V4_EtaPipPim_Lab, V4_EtaPip_Lab);
    Vector4<float>  V4_EtaPim_HeliEtaPipPim = helicityVec(V4_all_Lab, V4_EtaPipPim_Lab, V4_EtaPim_Lab);
    Vector4<float>  V4_PipPim_HeliEtaPipPim = helicityVec(V4_all_Lab, V4_EtaPipPim_Lab, V4_PipPim_Lab);
    
    Vector4<float>  V4_Pip_HeliPipPim = helicityVec(V4_EtaPipPim_Lab, V4_PipPim_Lab, V4_Pip_Lab);
    Vector4<float>  V4_Pim_HeliPipPim = helicityVec(V4_EtaPipPim_Lab, V4_PipPim_Lab, V4_Pim_Lab);
    
    Vector4<float>  V4_Pip_HeliEtaPip = helicityVec(V4_EtaPipPim_Lab, V4_EtaPip_Lab, V4_Pip_Lab);
    Vector4<float>  V4_Pim_HeliEtaPim = helicityVec(V4_EtaPipPim_Lab, V4_EtaPim_Lab, V4_Pim_Lab);
    Vector4<float>  V4_Eta_HeliEtaPip = helicityVec(V4_EtaPipPim_Lab, V4_EtaPip_Lab, V4_Eta_Lab);
    Vector4<float>  V4_Eta_HeliEtaPim = helicityVec(V4_EtaPipPim_Lab, V4_EtaPim_Lab, V4_Eta_Lab);
    
    Vector4<float> V4_normPipPimDecHeliEtaPipPim(0.5*(V4_Pip_HeliEtaPipPim.T()+V4_Pim_HeliEtaPipPim.T()),
						 V4_Pim_HeliEtaPipPim.Y()*V4_Pip_HeliEtaPipPim.Z()-V4_Pim_HeliEtaPipPim.Z()*V4_Pip_HeliEtaPipPim.Y(),
						 V4_Pim_HeliEtaPipPim.Z()*V4_Pip_HeliEtaPipPim.X()-V4_Pim_HeliEtaPipPim.X()*V4_Pip_HeliEtaPipPim.Z(),
						 V4_Pim_HeliEtaPipPim.X()*V4_Pip_HeliEtaPipPim.Y()-V4_Pim_HeliEtaPipPim.Y()*V4_Pip_HeliEtaPipPim.X());
    
    Vector4<float> V4_normEtaPipDecHeliEtaPipPim(0.5*(V4_Eta_HeliEtaPipPim.T()+V4_Pip_HeliEtaPipPim.T()),
						 V4_Pip_HeliEtaPipPim.Y()*V4_Eta_HeliEtaPipPim.Z()-V4_Pip_HeliEtaPipPim.Z()*V4_Eta_HeliEtaPipPim.Y(),
						 V4_Pip_HeliEtaPipPim.Z()*V4_Eta_HeliEtaPipPim.X()-V4_Pip_HeliEtaPipPim.X()*V4_Eta_HeliEtaPipPim.Z(),
						 V4_Pip_HeliEtaPipPim.X()*V4_Eta_HeliEtaPipPim.Y()-V4_Pip_HeliEtaPipPim.Y()*V4_Eta_HeliEtaPipPim.X());
    
    EvtDataNew* evtData=new EvtDataNew();
    evtData->FourVecsProd[enumProd4V::Psi] = V4_psi;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPipPim_HeliPsi] = V4_EtaPipPim_HeliPsi;
    
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::PipPim_HeliPsi] = V4_PipPim_HeliPsi;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPip_HeliPsi] = V4_EtaPip_HeliPsi;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPim_HeliPsi] = V4_EtaPim_HeliPsi;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::gamma_HeliPsi] = V4_gamma_HeliPsi;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::Eta_HeliPsi] = V4_Eta_HeliPsi;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::Pip_HeliPsi] = V4_Pip_HeliPsi;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::Pim_HeliPsi] = V4_Pim_HeliPsi;
    
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::Eta_HeliEtaPipPim] = V4_Eta_HeliEtaPipPim;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::Pip_HeliEtaPipPim] = V4_Pip_HeliEtaPipPim;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::Pim_HeliEtaPipPim] = V4_Pim_HeliEtaPipPim;
    
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPip_HeliEtaPipPim] = V4_EtaPip_HeliEtaPipPim;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPim_HeliEtaPipPim] = V4_EtaPim_HeliEtaPipPim;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::PipPim_HeliEtaPipPim] = V4_PipPim_HeliEtaPipPim;
    
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::Pip_HeliPipPim] = V4_Pip_HeliPipPim;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::Pim_HeliPipPim] = V4_Pim_HeliPipPim;
    
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::Pip_HeliEtaPip] = V4_Pip_HeliEtaPip;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::Pim_HeliEtaPim] = V4_Pim_HeliEtaPim;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::Eta_HeliEtaPip] = V4_Eta_HeliEtaPip;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::Eta_HeliEtaPim] = V4_Eta_HeliEtaPim;
    
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::normPipPimDecHeliEtaPipPim] = V4_normPipPimDecHeliEtaPipPim;
    evtData->FourVecsDec[enumJpsiGamEtaPiPi4V::normEtaPipDecHeliEtaPipPim] = V4_normEtaPipDecHeliEtaPipPim;
    
    // calculate and store WignerD functions for Psi -> X gamma (JPC =0-+, 1++, 2-+)

    Spin jPsi=1;
    for (Spin M=-1; M<=1; M=M+2){
      for (Spin lam=-1; lam<=1; lam++){
	evtData->WignerDsProd[enumProdDfunc::Psi][jPsi][M][lam]=Wigner_D(0.,V4_EtaPipPim_HeliPsi.Theta(),0,jPsi,M,lam);
      }
    }
    
    for (Spin J_etapipi=0; J_etapipi<=2; J_etapipi++){
      for (Spin lam_etapipi=-J_etapipi; lam_etapipi<=J_etapipi; lam_etapipi++){
	evtData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::etapipi][J_etapipi][lam_etapipi][0]
	  =Wigner_D(V4_normPipPimDecHeliEtaPipPim.Phi(),V4_normPipPimDecHeliEtaPipPim.Theta(),0,J_etapipi,lam_etapipi,0);
      }
    }
    
    for (Spin J_etapipi=0; J_etapipi<=2; J_etapipi++){
      for (Spin lam_etapipi=-J_etapipi; lam_etapipi<=J_etapipi; lam_etapipi++){
	for (Spin lamAorF=-J_etapipi; lamAorF<=J_etapipi; lamAorF++){
	  evtData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XToAplusPiminus][J_etapipi][lam_etapipi][lamAorF]
	    =Wigner_D(V4_EtaPip_HeliEtaPipPim.Phi(),V4_EtaPip_HeliEtaPipPim.Theta(),0,J_etapipi,lam_etapipi,lamAorF);
	  
	  evtData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XToAminusPiplus][J_etapipi][lam_etapipi][lamAorF]
	    =Wigner_D(V4_EtaPim_HeliEtaPipPim.Phi(),V4_EtaPim_HeliEtaPipPim.Theta(),0,J_etapipi,lam_etapipi,lamAorF);
	  
	  evtData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XTofEta][J_etapipi][lam_etapipi][lamAorF]
	    =Wigner_D(V4_PipPim_HeliEtaPipPim.Phi(), V4_PipPim_HeliEtaPipPim.Theta(),0,J_etapipi,lam_etapipi,lamAorF);
	  
	}
      }
    }
    
    for (Spin Ja=0; Ja<=2; Ja++){
      for (Spin lam_a=-Ja; lam_a<=Ja; lam_a++){
	evtData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::Aplus][Ja][lam_a][0]
	  =Wigner_D(V4_Pip_HeliEtaPip.Phi(), V4_Pip_HeliEtaPip.Theta(),0,Ja,lam_a,0);
	
	evtData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::Aminus][Ja][lam_a][0]
	  =Wigner_D(V4_Pim_HeliEtaPim.Phi(), V4_Pim_HeliEtaPim.Theta(),0,Ja,lam_a,0);
      }
    }
    
    for (Spin Jf=0; Jf<=2; Jf++){
      for (Spin lam_f=-Jf; lam_f<=Jf; lam_f++){
	evtData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::fToPiPi][Jf][lam_f][0]
	  =Wigner_D(V4_Pim_HeliPipPim.Phi(), V4_Pim_HeliPipPim.Theta(),0,Jf,lam_f,0);
	
      }
    }
    
    
    evtData->evtWeight=1.;
    theEvtList.push_back(evtData);
    
    ++evtCount;
  }
}


