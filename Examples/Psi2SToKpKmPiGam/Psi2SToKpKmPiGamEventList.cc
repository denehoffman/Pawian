#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEventList.hh"
#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


Psi2SToKpKmPiGamEventList::Psi2SToKpKmPiGamEventList(EventList& evtListData, EventList& evtListMc)
{
  read4Vecs(evtListData, _dataList);
  read4Vecs(evtListMc, _mcList);
}



Psi2SToKpKmPiGamEventList::~Psi2SToKpKmPiGamEventList()
{
}

void Psi2SToKpKmPiGamEventList::read4Vecs(EventList& evtList, std::vector<Psi2SToKpKmPiGamEvtData*>& theEvtList)
{

   Event* anEvent;
   int evtCount = 0;
   while ((anEvent = evtList.nextEvent()) != 0 && evtCount<30000){
//    while ((anEvent = evtList.nextEvent())){
     if ( evtCount%1000 == 0 ) Info << "4vec calculation for event " << evtCount << endmsg;
     Vector4<float> psi2S_4V(*(anEvent->p4(0))+*(anEvent->p4(1))+*(anEvent->p4(2))+*(anEvent->p4(3)));
     Vector4<float> chic1_4V(*(anEvent->p4(0))+*(anEvent->p4(1))+*(anEvent->p4(2)));     

     if ( evtCount%1000 == 0 ){
     Info << "Psi(2S) 4vec" << "\n"
	  << " px: " << psi2S_4V.Px() <<"\t"
	  << " py: " << psi2S_4V.Py() <<"\t"
	  << " pz: " << psi2S_4V.Pz() <<"\t" 
	  << " e : " << psi2S_4V.E() << "\t"
	  << " m : " << psi2S_4V.M() << endmsg;

     Info << "Chic1 4vec" << "\n"
	  << " px: " << chic1_4V.Px() <<"\t"
	  << " py: " << chic1_4V.Py() <<"\t"
	  << " pz: " << chic1_4V.Pz() <<"\t" 
	  << " e : " << chic1_4V.E() <<"\t"
 	  << " m : " << chic1_4V.M() << endmsg;
     }

     Vector4<float>  chic1_HeliPsi2S_4V(chic1_4V);
     chic1_HeliPsi2S_4V.Boost(psi2S_4V);   

     Vector4<float> KpKm_4V(*(anEvent->p4(0))+*(anEvent->p4(1)));
     Vector4<float> KpKm_HeliChic1_4V = helicityVec(psi2S_4V, chic1_4V, KpKm_4V);

     Vector4<float> KpPi_4V(*(anEvent->p4(0))+*(anEvent->p4(2)));
     Vector4<float> KpPi_HeliChic1_4V = helicityVec(psi2S_4V, chic1_4V, KpPi_4V);

     Vector4<float> KmPi_4V(*(anEvent->p4(1))+*(anEvent->p4(2)));
     Vector4<float> KmPi_HeliChic1_4V = helicityVec(psi2S_4V, chic1_4V, KmPi_4V);

     Vector4<float> Kp_4V(*(anEvent->p4(0)));
     Vector4<float> Kp_HeliKpKm_4V=helicityVec(chic1_4V, KpKm_4V, Kp_4V);
     Vector4<float> Kp_HeliKpPi_4V=helicityVec(chic1_4V, KpPi_4V, Kp_4V);

     Vector4<float> Km_4V(*(anEvent->p4(1)));
     Vector4<float> Km_HeliKmPi_4V=helicityVec(chic1_4V, KmPi_4V, Km_4V);

//        if (KpKm_4V.M()>3.12 && KpKm_4V.M()<3.22) continue;



     Psi2SToKpKmPiGamEvtData* thePsi2SToKpKmPiGamEvtData=new Psi2SToKpKmPiGamEvtData();
     thePsi2SToKpKmPiGamEvtData->psi2S_4V=psi2S_4V;
     thePsi2SToKpKmPiGamEvtData->chic1_HeliPsi2S_4V=chic1_HeliPsi2S_4V;
     thePsi2SToKpKmPiGamEvtData->KpKm_HeliChic1_4V=KpKm_HeliChic1_4V;
     thePsi2SToKpKmPiGamEvtData->KpPi_HeliChic1_4V=KpPi_HeliChic1_4V;
     thePsi2SToKpKmPiGamEvtData->KmPi_HeliChic1_4V=KmPi_HeliChic1_4V;
     thePsi2SToKpKmPiGamEvtData->Kp_HeliKpKm_4V=Kp_HeliKpKm_4V;
     thePsi2SToKpKmPiGamEvtData->Kp_HeliKpPi_4V=Kp_HeliKpPi_4V;
     thePsi2SToKpKmPiGamEvtData->Km_HeliKmPi_4V=Km_HeliKmPi_4V;

     // calculate and store WignerD functions for Psi(2S) -> Chi_c1 gamma
     Spin jPsi=1;
     for (Spin M=-1; M<=1; M=M+2){
 	for (Spin lam=-1; lam<=1; lam++){
	  thePsi2SToKpKmPiGamEvtData->DfPsi[jPsi][M][lam]=Wigner_D(chic1_HeliPsi2S_4V.Phi(),chic1_HeliPsi2S_4V.Theta(),0,jPsi,M,lam);
	}
      }

     // calculate and store WignerD functions for Chi_c1 -> K*+ K- , Chi_c1 -> K*- K+ and Chi_c1 -> a0(980) pi0
     Spin jChi_c1=1;
     for (Spin M=-1; M<=1; M++){
 	for (Spin lam=-1; lam<=1; lam++){
	  thePsi2SToKpKmPiGamEvtData->DfChiToKpPi[jChi_c1][M][lam]=Wigner_D(KpPi_HeliChic1_4V.Phi(),KpPi_HeliChic1_4V.Theta(),0,jChi_c1,M,lam);
	  thePsi2SToKpKmPiGamEvtData->DfChiToKmPi[jChi_c1][M][lam]=Wigner_D(KmPi_HeliChic1_4V.Phi(),KmPi_HeliChic1_4V.Theta(),0,jChi_c1,M,lam);
	  if (lam==0) thePsi2SToKpKmPiGamEvtData->DfChiToa0Pi[jChi_c1][M][lam]=Wigner_D(KpKm_HeliChic1_4V.Phi(),KpKm_HeliChic1_4V.Theta(),0,jChi_c1,M,lam);
	  thePsi2SToKpKmPiGamEvtData->DfChiToa2Pi[jChi_c1][M][lam]=Wigner_D(KpKm_HeliChic1_4V.Phi(),KpKm_HeliChic1_4V.Theta(),0,jChi_c1,M,lam);
	}
      }

     // calculate and store WignerD functions for K*1+ -> K+ pi0 and K*- -> K- pi0
     Spin jKstar1=1;
     for (Spin M=-1; M<=1; M++){
       Spin lam=0;
       thePsi2SToKpKmPiGamEvtData->DfKp1[jKstar1][M][lam]=Wigner_D(Kp_HeliKpPi_4V.Phi(),Kp_HeliKpPi_4V.Theta(),0,jKstar1,M,lam);
       thePsi2SToKpKmPiGamEvtData->DfKm1[jKstar1][M][lam]=Wigner_D(Km_HeliKmPi_4V.Phi(),Km_HeliKmPi_4V.Theta(),0,jKstar1,M,lam);
     }

     // calculate and store WignerD functions for K*2+ -> K+ pi0 and K*- -> K- pi0
     Spin jKstar2=2;
     for (Spin M=-2; M<=2; M++){
       Spin lam=0;
       thePsi2SToKpKmPiGamEvtData->DfKp2[jKstar2][M][lam]=Wigner_D(Kp_HeliKpPi_4V.Phi(),Kp_HeliKpPi_4V.Theta(),0,jKstar2,M,lam);
       thePsi2SToKpKmPiGamEvtData->DfKm2[jKstar2][M][lam]=Wigner_D(Km_HeliKmPi_4V.Phi(),Km_HeliKmPi_4V.Theta(),0,jKstar2,M,lam);
     }    

     // calculate and store WignerD functions for a2 -> K+ K- 
     Spin ja2=2;
     for (Spin M=-2; M<=2; M++){
       Spin lam=0;
       thePsi2SToKpKmPiGamEvtData->Dfa2[ja2][M][lam]=Wigner_D(Kp_HeliKpKm_4V.Phi(),Kp_HeliKpKm_4V.Theta(),0,ja2,M,lam);
     }

     theEvtList.push_back(thePsi2SToKpKmPiGamEvtData);

    ++evtCount;
   }
 
}


