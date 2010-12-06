#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


Psi2STo2K2PiGamEvtList::Psi2STo2K2PiGamEvtList(EventList& evtListData, EventList& evtListMc)
{ 
  Info << "read data events" << "\n";
  read4Vecs(evtListData, _dataList);

  Info << "read MC events" << "\n";
  read4Vecs(evtListMc, _mcList);
}



Psi2STo2K2PiGamEvtList::~Psi2STo2K2PiGamEvtList()
{
}

void Psi2STo2K2PiGamEvtList::read4Vecs(EventList& evtList, std::vector<Psi2STo2K2PiGamEvtData*>& theEvtList)
{

   Event* anEvent;
   int evtCount = 0;
   while ((anEvent = evtList.nextEvent()) != 0 && evtCount<50000){
//    while ((anEvent = evtList.nextEvent())){

     Vector4<float> psi2S_4V(*(anEvent->p4(0))+*(anEvent->p4(1))+*(anEvent->p4(2))+*(anEvent->p4(3))+*(anEvent->p4(4)));
     Vector4<float> chic0_4V(*(anEvent->p4(1))+*(anEvent->p4(2))+*(anEvent->p4(3))+*(anEvent->p4(4)));

     if ( evtCount%1000 == 0 ){
       Info << "4vec calculation for event " << evtCount << endmsg;
     
       Info << "Psi(2S) 4vec" << "\n"
	    << " px: " << psi2S_4V.Px() <<"\t"
	    << " py: " << psi2S_4V.Py() <<"\t"
	    << " pz: " << psi2S_4V.Pz() <<"\t" 
	    << " e : " << psi2S_4V.E() << "\t"
	    << " m : " << psi2S_4V.M() << endmsg;
       
       Info << "Chic0 4vec" << "\n"
	    << " px: " << chic0_4V.Px() <<"\t"
	    << " py: " << chic0_4V.Py() <<"\t"
	    << " pz: " << chic0_4V.Pz() <<"\t" 
	    << " e : " << chic0_4V.E() <<"\t"
	    << " m : " << chic0_4V.M() << endmsg;
     }

     Vector4<float>  chic0_HeliPsi2S_4V(chic0_4V);
     chic0_HeliPsi2S_4V.Boost(psi2S_4V);   

      Vector4<float> KpKm_4V(*(anEvent->p4(1))+*(anEvent->p4(2)));
      Vector4<float> KpKm_HeliChic0_4V = helicityVec(psi2S_4V, chic0_4V, KpKm_4V);

      Vector4<float> KpPi0_4V(*(anEvent->p4(1))+*(anEvent->p4(3)));
      Vector4<float> KpPi0_HeliChic0_4V = helicityVec(psi2S_4V, chic0_4V, KpPi0_4V);

      Vector4<float> KpPi1_4V(*(anEvent->p4(1))+*(anEvent->p4(4)));
      Vector4<float> KpPi1_HeliChic0_4V = helicityVec(psi2S_4V, chic0_4V, KpPi1_4V);

      Vector4<float> KmPi0_4V(*(anEvent->p4(2))+*(anEvent->p4(3)));
      Vector4<float> KmPi0_HeliChic0_4V = helicityVec(psi2S_4V, chic0_4V, KmPi0_4V);

      Vector4<float> KmPi1_4V(*(anEvent->p4(2))+*(anEvent->p4(4)));
      Vector4<float> KmPi1_HeliChic0_4V = helicityVec(psi2S_4V, chic0_4V, KmPi1_4V);

      Vector4<float> PiPi_4V(*(anEvent->p4(3))+*(anEvent->p4(4)));
      Vector4<float> PiPi_HeliChic0_4V = helicityVec(psi2S_4V, chic0_4V, PiPi_4V);

      Vector4<float> KpPiPi_4V(*(anEvent->p4(1))+*(anEvent->p4(3))+*(anEvent->p4(4)));
      Vector4<float> KpPiPi_HeliChic0_4V = helicityVec(psi2S_4V, chic0_4V, KpPiPi_4V);

      Vector4<float> KmPiPi_4V(*(anEvent->p4(2))+*(anEvent->p4(3))+*(anEvent->p4(4)));
      Vector4<float> KmPiPi_HeliChic0_4V = helicityVec(psi2S_4V, chic0_4V, KmPiPi_4V);

      Vector4<float> KKPi0_4V(*(anEvent->p4(1))+*(anEvent->p4(2))+*(anEvent->p4(3)));
      Vector4<float> KKPi0_HeliChic0_4V = helicityVec(psi2S_4V, chic0_4V, KKPi0_4V);

      Vector4<float> KKPi1_4V(*(anEvent->p4(1))+*(anEvent->p4(2))+*(anEvent->p4(4)));
      Vector4<float> KKPi1_HeliChic0_4V = helicityVec(psi2S_4V, chic0_4V, KKPi1_4V);

      Vector4<float> Kp_4V(*(anEvent->p4(1)));

//       Vector4<float> Kp_HeliKpPi0_4V=helicityVec(psi2S_4V, KpPi0_4V, Kp_4V);
//       Vector4<float> Kp_HeliKpPi1_4V=helicityVec(psi2S_4V, KpPi1_4V, Kp_4V);
//       Vector4<float> Kp_HeliKpPi0_ViaKpPiPi_4V=helicityVec(psi2S_4V, KpPi0_4V, Kp_4V);
//       Vector4<float> Kp_HeliKpPi1_ViaKpPiPi_4V=helicityVec(psi2S_4V, KpPi1_4V, Kp_4V);

      Vector4<float> Kp_HeliKpPi0_4V=helicityVec(chic0_4V, KpPi0_4V, Kp_4V);
      Vector4<float> Kp_HeliKpPi1_4V=helicityVec(chic0_4V, KpPi1_4V, Kp_4V);
      Vector4<float> Kp_HeliKpPi0_ViaKpPiPi_4V=helicityVec(KpPiPi_4V, KpPi0_4V, Kp_4V);
      Vector4<float> Kp_HeliKpPi1_ViaKpPiPi_4V=helicityVec(KpPiPi_4V, KpPi1_4V, Kp_4V);




      Vector4<float> Km_4V(*(anEvent->p4(2)));
      Vector4<float> Km_HeliKmPi0_4V=helicityVec(chic0_4V, KmPi0_4V, Km_4V);
      Vector4<float> Km_HeliKmPi1_4V=helicityVec(chic0_4V, KmPi1_4V, Km_4V);
      Vector4<float> Km_HeliKmPi0_ViaKmPiPi_4V=helicityVec(KmPiPi_4V, KmPi0_4V, Km_4V);
      Vector4<float> Km_HeliKmPi1_ViaKmPiPi_4V=helicityVec(KmPiPi_4V, KmPi1_4V, Km_4V);

//       Vector4<float> Km_HeliKmPi0_4V=helicityVec(psi2S_4V, KmPi0_4V, Km_4V);
//       Vector4<float> Km_HeliKmPi1_4V=helicityVec(psi2S_4V, KmPi1_4V, Km_4V);
//       Vector4<float> Km_HeliKmPi0_ViaKmPiPi_4V=helicityVec(psi2S_4V, KmPi0_4V, Km_4V);
//       Vector4<float> Km_HeliKmPi1_ViaKmPiPi_4V=helicityVec(psi2S_4V, KmPi1_4V, Km_4V);



      Vector4<float> KpPi0_HeliKpPi0Pi0_4V=helicityVec(chic0_4V, KpPiPi_4V, KpPi0_4V);
      Vector4<float> KpPi1_HeliKpPi0Pi0_4V=helicityVec(chic0_4V, KpPiPi_4V, KpPi1_4V);
      Vector4<float> KmPi0_HeliKmPi0Pi0_4V=helicityVec(chic0_4V, KmPiPi_4V, KmPi0_4V);
      Vector4<float> KmPi1_HeliKmPi0Pi0_4V=helicityVec(chic0_4V, KmPiPi_4V, KmPi1_4V);

//       Vector4<float> KpPi0_HeliKpPi0Pi0_4V=helicityVec(psi2S_4V, KpPiPi_4V, KpPi0_4V);
//       Vector4<float> KpPi1_HeliKpPi0Pi0_4V=helicityVec(psi2S_4V, KpPiPi_4V, KpPi1_4V);
//       Vector4<float> KmPi0_HeliKmPi0Pi0_4V=helicityVec(psi2S_4V, KmPiPi_4V, KmPi0_4V);
//       Vector4<float> KmPi1_HeliKmPi0Pi0_4V=helicityVec(psi2S_4V, KmPiPi_4V, KmPi1_4V);


      Vector4<float> PiPi_HeliKpPi0Pi0_4V=helicityVec(psi2S_4V, KpPiPi_4V, PiPi_4V);
      Vector4<float> PiPi_HeliKmPi0Pi0_4V=helicityVec(psi2S_4V, KmPiPi_4V, PiPi_4V);

      Vector4<float> Pi0_4V(*(anEvent->p4(3)));
      Vector4<float> Pi0_HeliPi0Pi0_ViaKpPi0Pi0_4V=helicityVec(KpPiPi_4V, PiPi_4V, Pi0_4V);
      Vector4<float> Pi0_HeliPi0Pi0_ViaKmPi0Pi0_4V=helicityVec(KmPiPi_4V, PiPi_4V, Pi0_4V);
//       Vector4<float> Pi0_HeliPi0Pi0_ViaKpPi0Pi0_4V=helicityVec(psi2S_4V, PiPi_4V, Pi0_4V);
//       Vector4<float> Pi0_HeliPi0Pi0_ViaKmPi0Pi0_4V=helicityVec(psi2S_4V, PiPi_4V, Pi0_4V);


      Psi2STo2K2PiGamEvtData* thePsi2STo2K2PiGamEvtData=new Psi2STo2K2PiGamEvtData();
      thePsi2STo2K2PiGamEvtData->psi2S_4V=psi2S_4V;
      thePsi2STo2K2PiGamEvtData->chic0_HeliPsi2S_4V=chic0_HeliPsi2S_4V;
      thePsi2STo2K2PiGamEvtData->KpKm_HeliChic0_4V=KpKm_HeliChic0_4V;
      thePsi2STo2K2PiGamEvtData->KpPi0_HeliChic0_4V=KpPi0_HeliChic0_4V;
      thePsi2STo2K2PiGamEvtData->KpPi1_HeliChic0_4V=KpPi1_HeliChic0_4V;
      thePsi2STo2K2PiGamEvtData->KmPi0_HeliChic0_4V=KmPi0_HeliChic0_4V;
      thePsi2STo2K2PiGamEvtData->KmPi1_HeliChic0_4V=KmPi1_HeliChic0_4V;
      thePsi2STo2K2PiGamEvtData->PiPi_HeliChic0_4V=PiPi_HeliChic0_4V;
      thePsi2STo2K2PiGamEvtData->KpPiPi_HeliChic0_4V=KpPiPi_HeliChic0_4V;
      thePsi2STo2K2PiGamEvtData->KmPiPi_HeliChic0_4V=KmPiPi_HeliChic0_4V;

      thePsi2STo2K2PiGamEvtData->KKPi0_HeliChic0_4V=KKPi0_HeliChic0_4V;
      thePsi2STo2K2PiGamEvtData->KKPi1_HeliChic0_4V=KKPi1_HeliChic0_4V;

      thePsi2STo2K2PiGamEvtData->KpPi0_HeliKpPi0Pi0_4V=KpPi0_HeliKpPi0Pi0_4V;
      thePsi2STo2K2PiGamEvtData->KpPi1_HeliKpPi0Pi0_4V=KpPi1_HeliKpPi0Pi0_4V;
      thePsi2STo2K2PiGamEvtData->KmPi0_HeliKmPi0Pi0_4V=KmPi0_HeliKmPi0Pi0_4V;
      thePsi2STo2K2PiGamEvtData->KmPi1_HeliKmPi0Pi0_4V=KmPi1_HeliKmPi0Pi0_4V;

      thePsi2STo2K2PiGamEvtData->Kp_HeliKpPi0_4V=Kp_HeliKpPi0_4V;
      thePsi2STo2K2PiGamEvtData->Kp_HeliKpPi1_4V=Kp_HeliKpPi1_4V;
      thePsi2STo2K2PiGamEvtData->Km_HeliKmPi0_4V=Km_HeliKmPi0_4V;
      thePsi2STo2K2PiGamEvtData->Km_HeliKmPi1_4V=Km_HeliKmPi1_4V;

      thePsi2STo2K2PiGamEvtData->Kp_HeliKpPi0_ViaKpPiPi_4V=Kp_HeliKpPi0_ViaKpPiPi_4V;
      thePsi2STo2K2PiGamEvtData->Kp_HeliKpPi1_ViaKpPiPi_4V=Kp_HeliKpPi1_ViaKpPiPi_4V;
      thePsi2STo2K2PiGamEvtData->Km_HeliKmPi0_ViaKmPiPi_4V=Km_HeliKmPi0_ViaKmPiPi_4V;
      thePsi2STo2K2PiGamEvtData->Km_HeliKmPi1_ViaKmPiPi_4V=Km_HeliKmPi1_ViaKmPiPi_4V;

      thePsi2STo2K2PiGamEvtData->PiPi_HeliKpPi0Pi0_4V=PiPi_HeliKpPi0Pi0_4V;
      thePsi2STo2K2PiGamEvtData->PiPi_HeliKmPi0Pi0_4V=PiPi_HeliKmPi0Pi0_4V;
      thePsi2STo2K2PiGamEvtData->Pi0_HeliPi0Pi0_ViaKpPi0Pi0_4V=Pi0_HeliPi0Pi0_ViaKpPi0Pi0_4V;
      thePsi2STo2K2PiGamEvtData->Pi0_HeliPi0Pi0_ViaKmPi0Pi0_4V=Pi0_HeliPi0Pi0_ViaKmPi0Pi0_4V;

     // calculate and store WignerD functions for Psi(2S) -> Chi_c0 gamma
     Spin jPsi=1;
     for (Spin M=-1; M<=1; M=M+2){
       for (Spin lam=-1; lam<=1; lam=lam+2){
	 thePsi2STo2K2PiGamEvtData->DfPsi[jPsi][M][lam]=Wigner_D(chic0_HeliPsi2S_4V.Phi(),chic0_HeliPsi2S_4V.Theta(),0,jPsi,M,lam);
       }
     }
     
 // calculate and store WignerD functions for K*+(1400)->K*+(892) pi1; K*+(892)->K+ pi0 + cc
     Spin jKst1400=1;
     for (Spin M=-1; M<=1; M++){
       for (Spin lam=-1; lam<=1; lam++){
	 thePsi2STo2K2PiGamEvtData->DfK1400pToKstpPi1[jKst1400][M][lam]=Wigner_D(KpPi0_HeliKpPi0Pi0_4V.Phi(),KpPi0_HeliKpPi0Pi0_4V.Theta(),0,jKst1400,M,lam);
	 thePsi2STo2K2PiGamEvtData->DfK1400pToKstpPi0[jKst1400][M][lam]=Wigner_D(KpPi1_HeliKpPi0Pi0_4V.Phi(),KpPi1_HeliKpPi0Pi0_4V.Theta(),0,jKst1400,M,lam);

	 thePsi2STo2K2PiGamEvtData->DfK1400mToKstmPi1[jKst1400][M][lam]=Wigner_D(KmPi0_HeliKmPi0Pi0_4V.Phi(),KmPi0_HeliKmPi0Pi0_4V.Theta(),0,jKst1400,M,lam);
	 thePsi2STo2K2PiGamEvtData->DfK1400mToKstmPi0[jKst1400][M][lam]=Wigner_D(KmPi1_HeliKmPi0Pi0_4V.Phi(),KmPi1_HeliKmPi0Pi0_4V.Theta(),0,jKst1400,M,lam);
       }
     }

 // calculate and store WignerD functions for K*1+->K+ pi0, K*1+->K+ pi1+ cc
     Spin jKst892=1;
     for (Spin M=-1; M<=1; M++){
       thePsi2STo2K2PiGamEvtData->DfKst1pToKpPi0[jKst892][M][0]=Wigner_D(Kp_HeliKpPi0_4V.Phi(),Kp_HeliKpPi0_4V.Theta(),0,jKst892,M,0);
       thePsi2STo2K2PiGamEvtData->DfKst1pToKpPi1[jKst892][M][0]=Wigner_D(Kp_HeliKpPi1_4V.Phi(),Kp_HeliKpPi1_4V.Theta(),0,jKst892,M,0);

       thePsi2STo2K2PiGamEvtData->DfKst1mToKmPi0[jKst892][M][0]=Wigner_D(Km_HeliKmPi0_4V.Phi(),Km_HeliKmPi0_4V.Theta(),0,jKst892,M,0);
       thePsi2STo2K2PiGamEvtData->DfKst1mToKmPi1[jKst892][M][0]=Wigner_D(Km_HeliKmPi1_4V.Phi(),Km_HeliKmPi1_4V.Theta(),0,jKst892,M,0);
     }

 // calculate and store WignerD functions for K*1+->K+ pi0, K*1+->K+ pi1+ cc Via K+ pi0 pi0
     for (Spin M=-1; M<=1; M++){
       thePsi2STo2K2PiGamEvtData->DfKst1pToKpPi0ViaKpPiPi[jKst892][M][0]=Wigner_D(Kp_HeliKpPi0_ViaKpPiPi_4V.Phi(), Kp_HeliKpPi0_ViaKpPiPi_4V.Theta(),0,jKst892,M,0);
       thePsi2STo2K2PiGamEvtData->DfKst1pToKpPi1ViaKpPiPi[jKst892][M][0]=Wigner_D(Kp_HeliKpPi1_ViaKpPiPi_4V.Phi(),Kp_HeliKpPi1_ViaKpPiPi_4V.Theta(),0,jKst892,M,0);

       thePsi2STo2K2PiGamEvtData->DfKst1mToKmPi0ViaKmPiPi[jKst892][M][0]=Wigner_D(Km_HeliKmPi0_ViaKmPiPi_4V.Phi(),Km_HeliKmPi0_ViaKmPiPi_4V.Theta(),0,jKst892,M,0);
       thePsi2STo2K2PiGamEvtData->DfKst1mToKmPi1ViaKmPiPi[jKst892][M][0]=Wigner_D(Km_HeliKmPi1_ViaKmPiPi_4V.Phi(),Km_HeliKmPi1_ViaKmPiPi_4V.Theta(),0,jKst892,M,0);
     }



 // calculate and store WignerD functions for K*2+->K+ pi0, K*2+->K+ pi1+ cc
     Spin jKst_2_1430=2;
     for (Spin M=-2; M<=2; M++){
       thePsi2STo2K2PiGamEvtData->DfKst2pToKpPi0[jKst_2_1430][M][0]=Wigner_D(Kp_HeliKpPi0_4V.Phi(),Kp_HeliKpPi0_4V.Theta(),0,jKst_2_1430,M,0);
       thePsi2STo2K2PiGamEvtData->DfKst2pToKpPi1[jKst_2_1430][M][0]=Wigner_D(Kp_HeliKpPi1_4V.Phi(),Kp_HeliKpPi1_4V.Theta(),0,jKst_2_1430,M,0);

       thePsi2STo2K2PiGamEvtData->DfKst2mToKmPi0[jKst_2_1430][M][0]=Wigner_D(Km_HeliKmPi0_4V.Phi(),Km_HeliKmPi0_4V.Theta(),0,jKst_2_1430,M,0);
       thePsi2STo2K2PiGamEvtData->DfKst2mToKmPi1[jKst_2_1430][M][0]=Wigner_D(Km_HeliKmPi1_4V.Phi(),Km_HeliKmPi1_4V.Theta(),0,jKst_2_1430,M,0);
     }

 

     theEvtList.push_back(thePsi2STo2K2PiGamEvtData);

    ++evtCount;
   }
 
}


