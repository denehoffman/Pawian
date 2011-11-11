#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEvtList.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


Psi2SToKpKmPiGamEvtList::Psi2SToKpKmPiGamEvtList(EventList& evtListData, EventList& evtListMc)
{
  read4Vecs(evtListData, _evtDataList);
  read4Vecs(evtListMc, _mcDataList);
}



Psi2SToKpKmPiGamEvtList::~Psi2SToKpKmPiGamEvtList()
{
}



void Psi2SToKpKmPiGamEvtList::read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList)
{

   Event* anEvent;
   int evtCount = 0;
   while ((anEvent = evtList.nextEvent()) != 0 && evtCount<30000){
//    while ((anEvent = evtList.nextEvent())){
     if ( evtCount%1000 == 0 ) Info << "4vec calculation for event " << evtCount << endmsg;
     Vector4<float> psi2S_4V(*(anEvent->p4(0))+*(anEvent->p4(1))+*(anEvent->p4(2))+*(anEvent->p4(3)));
     Vector4<double> psi2S_4Vd(&psi2S_4V);
     Vector4<float> chic1_4V(*(anEvent->p4(0))+*(anEvent->p4(1))+*(anEvent->p4(2)));     
     Vector4<double> chic1_4Vd(& chic1_4V);
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

     Vector4<double>  chic1_HeliPsi2S_4Vd(chic1_4Vd);
     chic1_HeliPsi2S_4Vd.Boost(psi2S_4Vd);

     Vector4<float> KpKm_4V(*(anEvent->p4(0))+*(anEvent->p4(1)));
     Vector4<double> KpKm_4Vd(&KpKm_4V);
     Vector4<double> KpKm_HeliChic1_4Vd = helicityVec(psi2S_4Vd, chic1_4Vd, KpKm_4Vd);


     Vector4<float> KpPi_4V(*(anEvent->p4(0))+*(anEvent->p4(2)));
     Vector4<double> KpPi_4Vd(&KpPi_4V);
     Vector4<double> KpPi_HeliChic1_4Vd = helicityVec(psi2S_4Vd, chic1_4Vd, KpPi_4Vd);

     Vector4<float> KmPi_4V(*(anEvent->p4(1))+*(anEvent->p4(2)));
     Vector4<double> KmPi_4Vd(&KmPi_4V);
     Vector4<double> KmPi_HeliChic1_4Vd = helicityVec(psi2S_4Vd, chic1_4Vd, KmPi_4Vd);

     Vector4<float> Kp_4V(*(anEvent->p4(0)));
     Vector4<double> Kp_4Vd(&Kp_4V);
     Vector4<double> Kp_HeliKpKm_4Vd=helicityVec(chic1_4Vd, KpKm_4Vd, Kp_4Vd);
     Vector4<double> Kp_HeliKpPi_4Vd=helicityVec(chic1_4Vd, KpPi_4Vd, Kp_4Vd);


     Vector4<float> Km_4V(*(anEvent->p4(1)));
     Vector4<double> Km_4Vd(&Km_4V);
     Vector4<double> Km_HeliKmPi_4Vd=helicityVec(chic1_4Vd, KmPi_4Vd, Km_4Vd);

     Vector4<float> Pi0_4V(*(anEvent->p4(2)));
     Vector4<double> Pi0_4Vd(&Pi0_4V);
     Vector4<double> Pi0_HeliChic1_4Vd = helicityVec(psi2S_4Vd, chic1_4Vd, Pi0_4Vd);
     Vector4<double> Kp_HeliChic1_4Vd = helicityVec(psi2S_4Vd, chic1_4Vd, Kp_4Vd);
     Vector4<double> Km_HeliChic1_4Vd = helicityVec(psi2S_4Vd, chic1_4Vd, Km_4Vd);

     Vector4<double> KpKmPi0_HeliChic1_4Vd(0.5*(Pi0_HeliChic1_4Vd.T()+Kp_HeliChic1_4Vd.T()+Km_HeliChic1_4Vd.T()),
                                            Pi0_HeliChic1_4Vd.Y()*Kp_HeliChic1_4Vd.Z()-Pi0_HeliChic1_4Vd.Z()*Kp_HeliChic1_4Vd.Y(),
                                            Pi0_HeliChic1_4Vd.Z()*Kp_HeliChic1_4Vd.X()-Pi0_HeliChic1_4Vd.X()*Kp_HeliChic1_4Vd.Z(),
                                            Pi0_HeliChic1_4Vd.X()*Kp_HeliChic1_4Vd.Y()-Pi0_HeliChic1_4Vd.Y()*Kp_HeliChic1_4Vd.X());


//        if (KpKm_4V.M()>3.12 && KpKm_4V.M()<3.22) continue;



     EvtData* evtData=new EvtData();
     evtData->FourVecs[enumChic1ToKpKmPiGamData::Psi2S]=psi2S_4V=psi2S_4Vd;
     evtData->FourVecs[enumChic1ToKpKmPiGamData::Chic1_HeliPsi2S]=chic1_HeliPsi2S_4Vd;
     evtData->FourVecs[enumChic1ToKpKmPiGamData::KpKm_HeliChic1]=KpKm_HeliChic1_4Vd;
     evtData->FourVecs[enumChic1ToKpKmPiGamData::KpPi_HeliChic1]=KpPi_HeliChic1_4Vd;
     evtData->FourVecs[enumChic1ToKpKmPiGamData::KmPi_HeliChic1]=KmPi_HeliChic1_4Vd;
     evtData->FourVecs[enumChic1ToKpKmPiGamData::KpKmPi0_HeliChic1]=KpKmPi0_HeliChic1_4Vd;
     evtData->FourVecs[enumChic1ToKpKmPiGamData::Kp_HeliKpKm]=Kp_HeliKpKm_4Vd;
     evtData->FourVecs[enumChic1ToKpKmPiGamData::Kp_HeliKpPi]=Kp_HeliKpPi_4Vd;
     evtData->FourVecs[enumChic1ToKpKmPiGamData::Km_HeliKmPi]=Km_HeliKmPi_4Vd;

     // calculate and store WignerD functions for Psi(2S) -> Chi_c1 gamma
     Spin jPsi=1;
     for (Spin M=-1; M<=1; M=M+2){
 	for (Spin lam=-1; lam<=1; lam++){
	  evtData->WignerDs[enumChic1ToKpKmPiGamData::dfPsi][jPsi][M][lam]=Wigner_D(0.,chic1_HeliPsi2S_4Vd.Theta(),0,jPsi,M,lam);
	}
     }

     // calculate and store WignerD functions for Chi_c1 -> K*+ K- , Chi_c1 -> K*- K+ and Chi_c1 -> a0(980) pi0
     Spin jChi_c1=1;
     for (Spin M=-1; M<=1; M++){
 	for (Spin lam=-1; lam<=1; lam++){
	  evtData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToKpPi][jChi_c1][M][lam]=Wigner_D(KpPi_HeliChic1_4Vd.Phi(),KpPi_HeliChic1_4Vd.Theta(),0,jChi_c1,M,lam);
	  evtData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToKmPi][jChi_c1][M][lam]=Wigner_D(KmPi_HeliChic1_4Vd.Phi(),KmPi_HeliChic1_4Vd.Theta(),0,jChi_c1,M,lam);
	  if (lam==0){
	    evtData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToa0Pi][jChi_c1][M][lam]=Wigner_D(KpKm_HeliChic1_4Vd.Phi(),KpKm_HeliChic1_4Vd.Theta(),0,jChi_c1,M,lam);
	    evtData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToKKPi][jChi_c1][M][lam]=Wigner_D(KpKmPi0_HeliChic1_4Vd.Phi(), KpKmPi0_HeliChic1_4Vd.Theta(),0,jChi_c1,M,lam);
	  }
	  evtData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToa2Pi][jChi_c1][M][lam]=Wigner_D(KpKm_HeliChic1_4Vd.Phi(),KpKm_HeliChic1_4Vd.Theta(),0,jChi_c1,M,lam);
	}
      }

     // calculate and store WignerD functions for K*1+ -> K+ pi0 and K*- -> K- pi0
     Spin jKstar1=1;
     for (Spin M=-1; M<=1; M++){
       Spin lam=0;
       evtData->WignerDs[enumChic1ToKpKmPiGamData::DfKp1][jKstar1][M][lam]=Wigner_D(Kp_HeliKpPi_4Vd.Phi(),Kp_HeliKpPi_4Vd.Theta(),0,jKstar1,M,lam);
       evtData->WignerDs[enumChic1ToKpKmPiGamData::DfKm1][jKstar1][M][lam]=Wigner_D(Km_HeliKmPi_4Vd.Phi(),Km_HeliKmPi_4Vd.Theta(),0,jKstar1,M,lam);
     }

     // calculate and store WignerD functions for K*2+ -> K+ pi0 and K*- -> K- pi0
     Spin jKstar2=2;
     for (Spin M=-2; M<=2; M++){
       Spin lam=0;
       evtData->WignerDs[enumChic1ToKpKmPiGamData::DfKp2][jKstar2][M][lam]=Wigner_D(Kp_HeliKpPi_4Vd.Phi(),Kp_HeliKpPi_4Vd.Theta(),0,jKstar2,M,lam);
       evtData->WignerDs[enumChic1ToKpKmPiGamData::DfKm2][jKstar2][M][lam]=Wigner_D(Km_HeliKmPi_4Vd.Phi(),Km_HeliKmPi_4Vd.Theta(),0,jKstar2,M,lam);
     }    

     // calculate and store WignerD functions for a2 -> K+ K- 
     Spin ja2=2;
     for (Spin M=-2; M<=2; M++){
       Spin lam=0;
       evtData->WignerDs[enumChic1ToKpKmPiGamData::Dfa2][ja2][M][lam]=Wigner_D(Kp_HeliKpKm_4Vd.Phi(),Kp_HeliKpKm_4Vd.Theta(),0,ja2,M,lam);
     }

     evtData->evtWeight=1.;
     theEvtList.push_back(evtData);

    ++evtCount;
   }
 
}


