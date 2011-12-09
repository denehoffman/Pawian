#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKEventList.hh"
#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


JpsiGamKsKlKKEventList::JpsiGamKsKlKKEventList(EventList& evtListData, EventList& evtListMc)
{
  read4Vecs(evtListData, _dataList);
  read4Vecs(evtListMc, _mcList);
}



JpsiGamKsKlKKEventList::~JpsiGamKsKlKKEventList()
{
}

void JpsiGamKsKlKKEventList::read4Vecs(EventList& evtList, std::vector<JpsiGamKsKlKKEvtData*>& theEvtList)
{

   Event* anEvent;
   int evtCount = 0;
   while ((anEvent = evtList.nextEvent())){
     if ( evtCount%1000 == 0 ) Info << "4vec calculation for event " << evtCount ;  // << endmsg;
     
     Vector4<float> gam = *(anEvent->p4(0));
     Vector4<float> ks = *(anEvent->p4(1));
     Vector4<float> kl = *(anEvent->p4(2));
     Vector4<float> kp = *(anEvent->p4(3));
     Vector4<float> km = *(anEvent->p4(4));
     
     Vector4<float> V4_psi = gam+ks+kl+kp+km;
    
     
     if ( evtCount%1000 == 0 ){
     Info << "psi 4vec" << "\n"
	  << " px: " << V4_psi.Px() <<"\t"
	  << " py: " << V4_psi.Py() <<"\t"
	  << " pz: " << V4_psi.Pz() <<"\t" 
	  << " e : " << V4_psi.E() << "\t"
	  << " m : " << V4_psi.M() ;  // << endmsg;
     }
     
     Vector4<float>  V4_KsKlKpKm_HeliPsi( ks+kl+kp+km   );
     V4_KsKlKpKm_HeliPsi.Boost(V4_psi);
     
     Vector4<float>  V4_gamma_HeliPsi( gam   );
     V4_gamma_HeliPsi.Boost(V4_psi);
     Vector4<float>  V4_KpKm_HeliPsi( kp+km  );
     V4_KpKm_HeliPsi.Boost(V4_psi);
     Vector4<float>  V4_KsKl_HeliPsi( ks+kl  );
     V4_KsKl_HeliPsi.Boost(V4_psi);
     
     Vector4<float>  V4_KsKl_HeliKsKlKpKm( ks+kl   );
     V4_KsKl_HeliKsKlKpKm.Boost( ks+kl+kp+km  );
     
     Vector4<float>  V4_KpKm_HeliKsKlKpKm( kp+km   );
     V4_KpKm_HeliKsKlKpKm.Boost( ks+kl+kp+km  );
     
     Vector4<float>  V4_Ks_HeliKsKl( ks );
     V4_Ks_HeliKsKl.Boost( ks+kl );
     
     Vector4<float>  V4_Kp_HeliKpKm( kp );
     V4_Kp_HeliKpKm.Boost( kp+km );
     
     

     JpsiGamKsKlKKEvtData* theJpsiGamKsKlKKEvtData=new JpsiGamKsKlKKEvtData();
     theJpsiGamKsKlKKEvtData->V4_Psi=V4_psi;
     theJpsiGamKsKlKKEvtData->V4_KsKlKpKm_HeliPsi = V4_KsKlKpKm_HeliPsi;
     
     theJpsiGamKsKlKKEvtData->V4_KsKlKpKm_HeliPsi = V4_KsKlKpKm_HeliPsi;
     theJpsiGamKsKlKKEvtData->V4_KsKl_HeliPsi = V4_KsKl_HeliPsi;
     theJpsiGamKsKlKKEvtData->V4_KpKm_HeliPsi = V4_KpKm_HeliPsi;
     theJpsiGamKsKlKKEvtData->V4_gamma_HeliPsi = V4_gamma_HeliPsi;
     
     theJpsiGamKsKlKKEvtData->V4_KsKl_HeliKsKlKpKm = V4_KsKl_HeliKsKlKpKm;
     theJpsiGamKsKlKKEvtData->V4_KpKm_HeliKsKlKpKm = V4_KpKm_HeliKsKlKpKm;
     theJpsiGamKsKlKKEvtData->V4_Ks_HeliKsKl = V4_Ks_HeliKsKl;
     theJpsiGamKsKlKKEvtData->V4_Kp_HeliKpKm = V4_Kp_HeliKpKm;


     


     // calculate and store WignerD functions for Psi -> X gamma (JPC =0-=, 0++, 2++)
     Spin jPsi=1;
     for (Spin M=-1; M<=1; M=M+2){
 	for (Spin lam=-1; lam<=1; lam++){
	  theJpsiGamKsKlKKEvtData->df_Psi[jPsi][M][lam]=Wigner_D(0.,V4_KsKlKpKm_HeliPsi.Theta(),0,jPsi,M,lam);
	}
     }
     
     //WignerD function for 2+ -> phi phi
     Spin jTensor =2;
     for(Spin M=-jTensor; M<=jTensor; M++){
       for (Spin lam=-2; lam<=2; lam++){
	 theJpsiGamKsKlKKEvtData->Df_Spin2[jTensor][M][lam] = Wigner_D(V4_KsKl_HeliKsKlKpKm.Phi(),V4_KsKl_HeliKsKlKpKm.Theta(), 0,jTensor,M,lam);
       }
     }
     
     //WignerD function for 0-/0+ -> phi phi
     Spin jScalar =0;
     Spin M =0;
     Spin lam =0;
     theJpsiGamKsKlKKEvtData->Df_Spin0[jScalar][M][lam] = Wigner_D(V4_KsKl_HeliKsKlKpKm.Phi(),V4_KsKl_HeliKsKlKpKm.Theta(), 0,jScalar,M,lam);
     
     
     //WignerD function for phi -> K+ K- and phi -> KS KL
     Spin phiSpin=1;
     for(Spin M=-phiSpin; M<=phiSpin; M++){
       Spin lam=0;
       theJpsiGamKsKlKKEvtData->Df_KsKl[phiSpin][M][lam] = Wigner_D(V4_Ks_HeliKsKl.Phi(),V4_Ks_HeliKsKl.Theta(), 0,phiSpin,M,lam);
       theJpsiGamKsKlKKEvtData->Df_KpKm[phiSpin][M][lam] = Wigner_D(V4_Kp_HeliKpKm.Phi(),V4_Kp_HeliKpKm.Theta(), 0,phiSpin,M,lam);
     }
     

     theEvtList.push_back(theJpsiGamKsKlKKEvtData);

    ++evtCount;
   }
 
}


