#include <getopt.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include "Examples/pbarp/pbarpEventList.hh"
#include "Examples/pbarp/pbarpEnv.hh"
#include "Examples/pbarp/IsobarDecay.hh"
#include "Examples/pbarp/IsobarDecayList.hh"

#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"
#include "Particle/Particle.hh"
#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Utils/PawianCollectionUtils.hh"


pbarpEventList::pbarpEventList() :
  EvtDataBaseListNew()
{
}


pbarpEventList::~pbarpEventList()
{
}

void pbarpEventList::read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList, double& evtWeightSum, int maxEvts){
  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = evtList.nextEvent())){
    if (evtCount>= maxEvts) break;
    if (evtCount%10000 == 0) Info << "4vec calculation for event " << evtCount ;  // << endmsg;

    Vector4<float> V4_all_lab(0.,0.,0.,0.);
    
    std::vector< Vector4<float> > finalState4Vecs;
    std::vector<Particle*>  finalStateParticles=pbarpEnv::instance()->finalStateParticles();
    std::map<std::string, Vector4<float> > particle4VecMap;
    std::vector<Particle*>::iterator itPart;
    int counter=0;
    for (itPart=finalStateParticles.begin(); itPart != finalStateParticles.end(); ++itPart){
      Vector4<float> current4Vec = *(anEvent->p4(counter));
      finalState4Vecs.push_back(current4Vec);
      particle4VecMap[ (*itPart)->name()]=current4Vec;
      V4_all_lab += current4Vec;
      counter++;
    }     
    
    if (evtCount%10000 == 0){
      Info << "4vec all in lab system" << "\n"
	   << " px: " << V4_all_lab.Px() <<"\t"
	   << " py: " << V4_all_lab.Py() <<"\t"
	   << " pz: " << V4_all_lab.Pz() <<"\t"
	   << " e : " << V4_all_lab.E() << "\t"
	   << " m : " << V4_all_lab.M() ;  // << endmsg;
    }

    std::vector<Particle*>  sortedFinalStateParticles=pbarpEnv::instance()->finalStateParticles();    
    pawian::Collection::PtrLess thePtrLess; 
    std::sort(sortedFinalStateParticles.begin(), sortedFinalStateParticles.end(), thePtrLess);    

    std::string name_all_lab=getName(sortedFinalStateParticles);


    EvtDataNew* evtData=new EvtDataNew();
    evtData->FourVecsString["all"]=V4_all_lab;
    //cache 4 vectors of inital state particles
    std::map<std::string, Vector4<float> >::iterator it4VecMap;
    for (it4VecMap=particle4VecMap.begin(); it4VecMap!=particle4VecMap.end(); ++it4VecMap){
      evtData->FourVecsString[it4VecMap->first]=it4VecMap->second;
    }
 
   //fill WignerD functions
    std::vector<boost::shared_ptr<IsobarDecay> > theDecays=pbarpEnv::instance()->prodDecayList()->getList();
    std::vector<boost::shared_ptr<IsobarDecay> >::iterator itIso;
    for (itIso=theDecays.begin(); itIso!=theDecays.end(); ++itIso){
      (*itIso)->fillWignerDs(particle4VecMap, evtData);
    }
 
    // Vector4<float>  V4_all_Lab( ks+kl+kp+km+gam   );
    //  Vector4<float>  V4_KsKlKpKm_Lab( ks+kl+kp+km   );
    //  Vector4<float>  V4_KsKl_Lab( ks+kl   );
    //  Vector4<float>  V4_Ks_Lab( ks   );
    //  Vector4<float>  V4_KpKm_Lab( kp+km   );
    //  Vector4<float>  V4_Kp_Lab( kp   );
    //  Vector4<float>  V4_Km_Lab( km   );
    //  Vector4<float>  V4_Kl_Lab( kl   );
          
    //  Vector4<float>  V4_KsKlKpKm_HeliPsi( ks+kl+kp+km   );
    //  V4_KsKlKpKm_HeliPsi.Boost(V4_psi);
     
    //  Vector4<float>  V4_gamma_HeliPsi( gam   );
    //  V4_gamma_HeliPsi.Boost(V4_psi);
    //  Vector4<float>  V4_KpKm_HeliPsi( kp+km  );
    //  V4_KpKm_HeliPsi.Boost(V4_psi);
    //  Vector4<float>  V4_KsKl_HeliPsi( ks+kl  );
    //  V4_KsKl_HeliPsi.Boost(V4_psi);

    //  Vector4<float>  V4_Kl_HeliPsi( kl   );
    //  V4_Kl_HeliPsi.Boost( V4_psi );
    //  Vector4<float>  V4_Ks_HeliPsi( ks   );
    //  V4_Ks_HeliPsi.Boost( V4_psi );
    //  Vector4<float>  V4_Kp_HeliPsi( kp   );
    //  V4_Kp_HeliPsi.Boost( V4_psi );
    //  Vector4<float>  V4_Km_HeliPsi( km   );
    //  V4_Km_HeliPsi.Boost( V4_psi );

    //  Vector4<float>  V4_KsKl_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_KsKl_Lab);
    //  Vector4<float>  V4_KpKm_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_KpKm_Lab);     
    //  Vector4<float>  V4_Ks_HeliKsKl=helicityVec(V4_KsKlKpKm_Lab, V4_KsKl_Lab, V4_Ks_Lab);     
    //  Vector4<float>  V4_Kp_HeliKpKm=helicityVec(V4_KsKlKpKm_Lab, V4_KpKm_Lab, V4_Kp_Lab);

    //  Vector4<float>  V4_Kp_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_Kp_Lab);
    //  Vector4<float>  V4_Km_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_Km_Lab);
     
    //  Vector4<float> V4_normKpKmDecHeliKsKlKpKm
    //    (0.5*(V4_Kp_HeliKsKlKpKm.T()+V4_Km_HeliKsKlKpKm.T()),
    // 	V4_Km_HeliKsKlKpKm.Y()*V4_Kp_HeliKsKlKpKm.Z()-V4_Km_HeliKsKlKpKm.Z()*V4_Kp_HeliKsKlKpKm.Y(),
    // 	V4_Km_HeliKsKlKpKm.Z()*V4_Kp_HeliKsKlKpKm.X()-V4_Km_HeliKsKlKpKm.X()*V4_Kp_HeliKsKlKpKm.Z(),
    // 	V4_Km_HeliKsKlKpKm.X()*V4_Kp_HeliKsKlKpKm.Y()-V4_Km_HeliKsKlKpKm.Y()*V4_Kp_HeliKsKlKpKm.X());
     
    //  Vector4<float>  V4_Kl_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_Kl_Lab);
    //  Vector4<float>  V4_Ks_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_Ks_Lab);
     
    //  Vector4<float> V4_normKsKlDecHeliKsKlKpKm
    //    (0.5*(V4_Kl_HeliKsKlKpKm.T()+V4_Ks_HeliKsKlKpKm.T()),
    // 	V4_Ks_HeliKsKlKpKm.Y()*V4_Kl_HeliKsKlKpKm.Z()-V4_Ks_HeliKsKlKpKm.Z()*V4_Kl_HeliKsKlKpKm.Y(),
    // 	V4_Ks_HeliKsKlKpKm.Z()*V4_Kl_HeliKsKlKpKm.X()-V4_Ks_HeliKsKlKpKm.X()*V4_Kl_HeliKsKlKpKm.Z(),
    // 	V4_Ks_HeliKsKlKpKm.X()*V4_Kl_HeliKsKlKpKm.Y()-V4_Ks_HeliKsKlKpKm.Y()*V4_Kl_HeliKsKlKpKm.X());
     
 
    //     EvtDataNew* evtData=new EvtDataNew();
     evtData->FourVecsString[name_all_lab] = V4_all_lab;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_KsKlKpKm_HeliPsi]=V4_KsKlKpKm_HeliPsi;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_KsKl_HeliPsi]=V4_KsKl_HeliPsi;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_KpKm_HeliPsi]=V4_KpKm_HeliPsi;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_gamma_HeliPsi]=V4_gamma_HeliPsi;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_Ks_HeliPsi]=V4_Ks_HeliPsi;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_Kl_HeliPsi]=V4_Kl_HeliPsi;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_Kp_HeliPsi]=V4_Kp_HeliPsi;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_Km_HeliPsi]=V4_Km_HeliPsi;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_KsKl_HeliKsKlKpKm]=V4_KsKl_HeliKsKlKpKm;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_KpKm_HeliKsKlKpKm]=V4_KpKm_HeliKsKlKpKm;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_Ks_HeliKsKl]=V4_Ks_HeliKsKl;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_Kp_HeliKpKm]= V4_Kp_HeliKpKm;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_normKpKmDecHeliKsKlKpKm]=V4_normKpKmDecHeliKsKlKpKm;
     // evtData->FourVecsDec[enumJpsiGamX4V::V4_normKsKlDecHeliKsKlKpKm]=V4_normKsKlDecHeliKsKlKpKm;

    
  //   // calculate and store WignerD functions for Psi -> X gamma (JPC =0-+, 0++, 2++)

  //   Spin jPsi=1;
  //   for (Spin M=-1; M<=1; M=M+2){
  //     for (Spin lam=-1; lam<=1; lam++){
  // 	evtData->WignerDsProd[enumProdDfunc::Psi][jPsi][M][lam]=Wigner_D(0.,V4_KsKlKpKm_HeliPsi.Theta(),0,jPsi,M,lam);
  //     }
  //   }

  // //WignerD functions for X -> phi phi

  //   for (Spin J_X=0; J_X<=2; J_X++){
  //     for (Spin lam_X=-J_X; lam_X<=J_X; lam_X++){
  // 	for (Spin lamPhi1mlamPhi2=-J_X; lamPhi1mlamPhi2<=J_X; lamPhi1mlamPhi2++){
  // 	evtData->WignerDsDec[enumJpsiGamXDfunc::Df_XToPhiPhi1][J_X][lam_X][lamPhi1mlamPhi2]
  // 	  =Wigner_D(V4_KsKl_HeliKsKlKpKm.Phi(),V4_KsKl_HeliKsKlKpKm.Theta(),0,J_X,lam_X,lamPhi1mlamPhi2);
  // 	evtData->WignerDsDec[enumJpsiGamXDfunc::Df_XToPhiPhi2][J_X][lam_X][lamPhi1mlamPhi2]
  // 	  =Wigner_D(V4_KpKm_HeliKsKlKpKm.Phi(),V4_KpKm_HeliKsKlKpKm.Theta(),0,J_X,lam_X,lamPhi1mlamPhi2);
  //     }
  //   }
  //   }
    
  //   //WignerD function for phi -> K+ K- and phi -> KS KL
  //   Spin phiSpin=1;
  //   for(Spin M=-phiSpin; M<=phiSpin; M++){
  //     Spin lam=0;
  //     evtData->WignerDsDec[enumJpsiGamXDfunc::Df_PhiToKsKl][phiSpin][M][lam] = Wigner_D(V4_Ks_HeliKsKl.Phi(),V4_Ks_HeliKsKl.Theta(), 0,phiSpin,M,lam);
  //     evtData->WignerDsDec[enumJpsiGamXDfunc::Df_PhiToKpKm][phiSpin][M][lam] = Wigner_D(V4_Kp_HeliKpKm.Phi(),V4_Kp_HeliKpKm.Theta(), 0,phiSpin,M,lam);
  //    }  
    
    
    evtData->evtWeight=anEvent->Weight();
    evtData->evtNo=_evtNoAll;
    theEvtList.push_back(evtData);
    
    evtWeightSum += anEvent->Weight();    
    ++evtCount;
    ++_evtNoAll;
  }
}

std::string pbarpEventList::getName(std::vector<Particle*>& theVec){
  std::string result;

  std::vector<Particle*>::iterator it;
  for(it=theVec.begin(); it!=theVec.end(); ++it){
    result+=(*it)->name();
  }

  return result;
}
