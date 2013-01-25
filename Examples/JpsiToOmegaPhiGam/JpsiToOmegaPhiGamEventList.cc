#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiToOmegaPhiGam/JpsiToOmegaPhiGamEventList.hh"

#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


JpsiToOmegaPhiGamEventList::JpsiToOmegaPhiGamEventList() :
  EvtDataBaseListNew()
{
}


JpsiToOmegaPhiGamEventList::~JpsiToOmegaPhiGamEventList()
{
}

void JpsiToOmegaPhiGamEventList::read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList, double& evtWeightSum, int maxEvts){
  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = evtList.nextEvent())){
    if (evtCount>= maxEvts) break;
    if (evtCount%10000 == 0) Info << "4vec calculation for event " << evtCount ;  // << endmsg;

    Vector4<float> gam = *(anEvent->p4(0));
    Vector4<float> kp = *(anEvent->p4(1));
    Vector4<float> km = *(anEvent->p4(2));
    Vector4<float> pip = *(anEvent->p4(3));
    Vector4<float> pim = *(anEvent->p4(4));
    Vector4<float> pi0 = *(anEvent->p4(5));
    
    Vector4<float> V4_psi = gam+kp+km+pip+pim+pi0;
    
    if (evtCount%10000 == 0){
      Info << "psi 4vec" << "\n"
	   << " px: " << V4_psi.Px() <<"\t"
	   << " py: " << V4_psi.Py() <<"\t"
	   << " pz: " << V4_psi.Pz() <<"\t"
	   << " e : " << V4_psi.E() << "\t"
	   << " m : " << V4_psi.M() ;  // << endmsg;
    }
    
     Vector4<float>  V4_all_Lab( pip+pim+pi0+kp+km+gam   );
     Vector4<float>  V4_PipPimPi0KpKm_Lab( pip+pim+pi0+kp+km   );
     Vector4<float>  V4_PipPimPi0_Lab( pip+pim+pi0   );
     Vector4<float>  V4_pip_Lab( pip   );
     Vector4<float>  V4_KpKm_Lab( kp+km   );
     Vector4<float>  V4_Kp_Lab( kp   );
     Vector4<float>  V4_Km_Lab( km   );
     Vector4<float>  V4_pim_Lab( pim   );
     Vector4<float>  V4_pi0_Lab( pi0   );
     
     Vector4<float>  V4_PipPimPi0KpKm_HeliPsi( V4_PipPimPi0KpKm_Lab   );
     V4_PipPimPi0KpKm_HeliPsi.Boost(V4_psi);
     
     Vector4<float>  V4_gamma_HeliPsi( gam   );
     V4_gamma_HeliPsi.Boost(V4_psi);
     Vector4<float>  V4_KpKm_HeliPsi( kp+km  );
     V4_KpKm_HeliPsi.Boost(V4_psi);
     Vector4<float>  V4_PipPimPi0_HeliPsi( pip+pim+pi0  );
     V4_PipPimPi0_HeliPsi.Boost(V4_psi);

     Vector4<float>  V4_Pip_HeliPsi( pip   );
     V4_Pip_HeliPsi.Boost( V4_psi );
     Vector4<float>  V4_Pim_HeliPsi( pim   );
     V4_Pim_HeliPsi.Boost( V4_psi );
     Vector4<float>  V4_Pi0_HeliPsi( pi0   );
     V4_Pi0_HeliPsi.Boost( V4_psi );
     Vector4<float>  V4_Kp_HeliPsi( kp   );
     V4_Kp_HeliPsi.Boost( V4_psi );
     Vector4<float>  V4_Km_HeliPsi( km   );
     V4_Km_HeliPsi.Boost( V4_psi );
     
     Vector4<float>  V4_PipPimPi0_HeliPipPimPi0KpKm=helicityVec(V4_all_Lab, V4_PipPimPi0KpKm_Lab, V4_PipPimPi0_Lab);
     Vector4<float>  V4_KpKm_HeliPipPimPi0KpKm=helicityVec(V4_all_Lab, V4_PipPimPi0KpKm_Lab, V4_KpKm_Lab);     
     Vector4<float>  V4_Kp_HeliPipPimPi0KpKm=helicityVec(V4_all_Lab, V4_PipPimPi0KpKm_Lab, V4_Kp_Lab);
     Vector4<float>  V4_Km_HeliPipPimPi0KpKm=helicityVec(V4_all_Lab, V4_PipPimPi0KpKm_Lab, V4_Km_Lab);
     
     Vector4<float>  V4_Pip_HeliPipPimPi0=helicityVec(V4_PipPimPi0KpKm_Lab, V4_PipPimPi0_Lab, V4_pip_Lab);
     Vector4<float>  V4_Pim_HeliPipPimPi0=helicityVec(V4_PipPimPi0KpKm_Lab, V4_PipPimPi0_Lab, V4_pim_Lab);
     Vector4<float>  V4_Pi0_HeliPipPimPi0=helicityVec(V4_PipPimPi0KpKm_Lab, V4_PipPimPi0_Lab, V4_pi0_Lab);
     
     Vector4<float>  V4_Kp_HeliKpKm=helicityVec(V4_PipPimPi0KpKm_Lab, V4_KpKm_Lab, V4_Kp_Lab);
     
     
     //
     // Omega decay plane normal and lambda parameter
     //
     Vector4<float> V4_omegaDecPlaneNormal_HeliPipPimPi0 ( 0.5*( V4_Pip_HeliPipPimPi0.T()+ V4_Pim_HeliPipPimPi0.T()+ V4_Pi0_HeliPipPimPi0.T() ),
							     V4_Pim_HeliPipPimPi0.Y()*V4_Pip_HeliPipPimPi0.Z() - V4_Pim_HeliPipPimPi0.Z()*V4_Pip_HeliPipPimPi0.Y(),
							     V4_Pim_HeliPipPimPi0.Z()*V4_Pip_HeliPipPimPi0.X() - V4_Pim_HeliPipPimPi0.X()*V4_Pip_HeliPipPimPi0.Z(),
							     V4_Pim_HeliPipPimPi0.X()*V4_Pip_HeliPipPimPi0.Y() - V4_Pim_HeliPipPimPi0.Y()*V4_Pip_HeliPipPimPi0.X() );
     double theQ =  V4_Pip_HeliPipPimPi0.E()- V4_Pip_HeliPipPimPi0.M() +  V4_Pim_HeliPipPimPi0.E()- V4_Pim_HeliPipPimPi0.M() +  V4_Pi0_HeliPipPimPi0.E()- V4_Pi0_HeliPipPimPi0.M();
     double lambdaNorm=theQ*theQ*(theQ*theQ/108.+ V4_Pim_HeliPipPimPi0.M()*theQ/9.+ V4_Pim_HeliPipPimPi0.M()* V4_Pim_HeliPipPimPi0.M()/3.);
     double lambda= V4_omegaDecPlaneNormal_HeliPipPimPi0.P()*  V4_omegaDecPlaneNormal_HeliPipPimPi0.P()/lambdaNorm;
     ///////
     


     //
     // omega and phi decay plane normal in omega phi helicity system
     //
     
     //omega
     Vector4<float> normalpipluspiminusheliomegaTLVlab= helicityVecInverse( V4_PipPimPi0KpKm_Lab, V4_PipPimPi0_Lab, V4_omegaDecPlaneNormal_HeliPipPimPi0);
     Vector4<float> V4_omegaDecPlaneNormal_HeliPipPimPi0KpKm= helicityVec(V4_all_Lab,V4_PipPimPi0KpKm_Lab , normalpipluspiminusheliomegaTLVlab);
     
     //test
     Vector4<float> testBoosted =  helicityVec(  V4_PipPimPi0KpKm_Lab, V4_PipPimPi0_Lab, normalpipluspiminusheliomegaTLVlab );
     //cout << "++++++++++++++++++++++++++++++"<<endl;
     //cout << V4_omegaDecPlaneNormal_HeliPipPimPi0 << endl;
     //cout <<  testBoosted  << endl;
     //cout << "+++++++++++++++++++++++++++++++" << endl;
     //eot
     



     
     //phi
     Vector4<float> V4_phiDecPlaneNormal_HeliPipPimPi0KpKm ( 0.5*( V4_Kp_HeliPipPimPi0KpKm.T()+ V4_Km_HeliPipPimPi0KpKm.T() ),
							V4_Km_HeliPipPimPi0KpKm.Y()*V4_Kp_HeliPipPimPi0KpKm.Z() - V4_Km_HeliPipPimPi0KpKm.Z()*V4_Kp_HeliPipPimPi0KpKm.Y(),
							V4_Km_HeliPipPimPi0KpKm.Z()*V4_Kp_HeliPipPimPi0KpKm.X() - V4_Km_HeliPipPimPi0KpKm.X()*V4_Kp_HeliPipPimPi0KpKm.Z(),
							V4_Km_HeliPipPimPi0KpKm.X()*V4_Kp_HeliPipPimPi0KpKm.Y() - V4_Km_HeliPipPimPi0KpKm.Y()*V4_Kp_HeliPipPimPi0KpKm.X() );
     //
     
     EvtDataNew* evtData=new EvtDataNew();
     evtData->FourVecsProd[enumProd4V::Psi] = V4_psi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_PipPimPi0KpKm_HeliPsi]=V4_PipPimPi0KpKm_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_PipPimPi0_HeliPsi]=V4_PipPimPi0_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_KpKm_HeliPsi]=V4_KpKm_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_gamma_HeliPsi]=V4_gamma_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Pip_HeliPsi]=V4_Pip_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Pim_HeliPsi]=V4_Pim_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Pi0_HeliPsi]=V4_Pi0_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Kp_HeliPsi]=V4_Kp_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Km_HeliPsi]=V4_Km_HeliPsi;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_PipPimPi0_HeliPipPimPi0KpKm]=V4_PipPimPi0_HeliPipPimPi0KpKm;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_KpKm_HeliPipPimPi0KpKm]=V4_KpKm_HeliPipPimPi0KpKm;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Kp_HeliKpKm]= V4_Kp_HeliKpKm;
     
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Pip_HeliPipPimPi0]=V4_Pip_HeliPipPimPi0;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Pim_HeliPipPimPi0]=V4_Pim_HeliPipPimPi0;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_Pi0_HeliPipPimPi0]=V4_Pi0_HeliPipPimPi0;
     
     evtData->FourVecsDec[enumJpsiGamX4V::V4_omegaDecPlaneNormal_HeliPipPimPi0] =V4_omegaDecPlaneNormal_HeliPipPimPi0;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_omegaDecPlaneNormal_HeliPipPimPi0KpKm] =V4_omegaDecPlaneNormal_HeliPipPimPi0KpKm;
     evtData->FourVecsDec[enumJpsiGamX4V::V4_phiDecPlaneNormal_HeliPipPimPi0KpKm] =V4_phiDecPlaneNormal_HeliPipPimPi0KpKm;
     
     evtData->KinematicVariables[enumJpsiGamXKin::OmegaDecLambda]=lambda;
     
     //
     // calculate and store WignerD functions for Psi -> X gamma (JPC =0-+, 0++, 2++)
     //
    Spin jPsi=1;
    for (Spin M=-1; M<=1; M=M+2){
      for (Spin lam=-1; lam<=1; lam++){
	evtData->WignerDsProd[enumProdDfunc::Psi][jPsi][M][lam]=Wigner_D(0.,V4_PipPimPi0KpKm_HeliPsi.Theta(),0,jPsi,M,lam);
      }
    }
    
  // //WignerD functions for X -> phi phi
    
    for (Spin J_X=0; J_X<=2; J_X++){
      for (Spin lam_X=-J_X; lam_X<=J_X; lam_X++){
   	for (Spin lamPhilamOmega=-J_X; lamPhilamOmega<=J_X; lamPhilamOmega++){
	  evtData->WignerDsDec[enumJpsiGamXDfunc::Df_XToOmegaPhi_KK][J_X][lam_X][lamPhilamOmega]
	    =Wigner_D(V4_KpKm_HeliPipPimPi0KpKm.Phi(),V4_KpKm_HeliPipPimPi0KpKm.Theta(),0,J_X,lam_X,lamPhilamOmega);
	  evtData->WignerDsDec[enumJpsiGamXDfunc::Df_XToOmegaPhi_3pi][J_X][lam_X][lamPhilamOmega]
	    =Wigner_D(V4_PipPimPi0_HeliPipPimPi0KpKm.Phi(), V4_PipPimPi0_HeliPipPimPi0KpKm.Theta(),0,J_X,lam_X,lamPhilamOmega);
	}
      }
    }
    
    //WignerD function for phi -> K+ K- andomega -> pi+ pi- pi0
    Spin phiSpin=1;
    for(Spin M=-phiSpin; M<=phiSpin; M++){
      Spin lam=0;
      evtData->WignerDsDec[enumJpsiGamXDfunc::Df_PhiToKpKm][phiSpin][M][lam] = Wigner_D(V4_Kp_HeliKpKm.Phi(),V4_Kp_HeliKpKm.Theta(), 0,phiSpin,M,lam);
      evtData->WignerDsDec[enumJpsiGamXDfunc::Df_OmegaToPipPimPi0][phiSpin][M][lam] = Wigner_D(V4_omegaDecPlaneNormal_HeliPipPimPi0.Phi(), V4_omegaDecPlaneNormal_HeliPipPimPi0.Theta(), 0,phiSpin,M,lam);
    }    
    
    evtData->evtWeight=anEvent->Weight();
    evtData->evtNo=_evtNoAll;
    theEvtList.push_back(evtData);
    
    evtWeightSum += anEvent->Weight();    
    ++evtCount;
    ++_evtNoAll;
  }
}


