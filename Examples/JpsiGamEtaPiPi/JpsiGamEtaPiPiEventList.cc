#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiEventList.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiData.hh"

#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


JpsiGamEtaPiPiEventList::JpsiGamEtaPiPiEventList(EventList& evtListData, EventList& evtListMc)
{
  read4Vecs(evtListData, _evtDataList);
  read4Vecs(evtListMc, _mcDataList);
}


JpsiGamEtaPiPiEventList::~JpsiGamEtaPiPiEventList()
{
}

void JpsiGamEtaPiPiEventList::read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList){
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

		Vector4<float>  V4_Pip_HeliPipPim = helicityVec(V4_all_Lab, V4_PipPim_Lab, V4_Pip_Lab);
		Vector4<float>  V4_Pim_HeliPipPim = helicityVec(V4_all_Lab, V4_PipPim_Lab, V4_Pim_Lab);

		Vector4<float>  V4_Pip_HeliEtaPip = helicityVec(V4_all_Lab, V4_EtaPip_Lab, V4_Pip_Lab);
		Vector4<float>  V4_Pim_HeliEtaPim = helicityVec(V4_all_Lab, V4_EtaPim_Lab, V4_Pim_Lab);
		Vector4<float>  V4_Eta_HeliEtaPip = helicityVec(V4_all_Lab, V4_EtaPip_Lab, V4_Eta_Lab);
		Vector4<float>  V4_Eta_HeliEtaPim = helicityVec(V4_all_Lab, V4_EtaPim_Lab, V4_Eta_Lab);
/*
		Vector4<float>  V4_EtaPip_HeliEtaPipPim=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_KsKl_Lab);
		Vector4<float>  V4_KpKm_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_KpKm_Lab);
		Vector4<float>  V4_Ks_HeliKsKl=helicityVec(V4_KsKlKpKm_Lab, V4_KsKl_Lab, V4_Ks_Lab);
		Vector4<float>  V4_Kp_HeliKpKm=helicityVec(V4_KsKlKpKm_Lab, V4_KpKm_Lab, V4_Kp_Lab);

		Vector4<float>  V4_Kp_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_Kp_Lab);
		Vector4<float>  V4_Km_HeliKsKlKpKm=helicityVec(V4_all_Lab, V4_KsKlKpKm_Lab, V4_Km_Lab);
*/
		Vector4<float> V4_normPipPimDecHeliEtaPipPim(0.5*(V4_Pip_HeliEtaPipPim.T()+V4_Pim_HeliEtaPipPim.T()),
			V4_Pim_HeliEtaPipPim.Y()*V4_Pip_HeliEtaPipPim.Z()-V4_Pim_HeliEtaPipPim.Z()*V4_Pip_HeliEtaPipPim.Y(),
			V4_Pim_HeliEtaPipPim.Z()*V4_Pip_HeliEtaPipPim.X()-V4_Pim_HeliEtaPipPim.X()*V4_Pip_HeliEtaPipPim.Z(),
			V4_Pim_HeliEtaPipPim.X()*V4_Pip_HeliEtaPipPim.Y()-V4_Pim_HeliEtaPipPim.Y()*V4_Pip_HeliEtaPipPim.X());

		Vector4<float> V4_normEtaPipDecHeliEtaPipPim(0.5*(V4_Eta_HeliEtaPipPim.T()+V4_Pip_HeliEtaPipPim.T()),
			V4_Pip_HeliEtaPipPim.Y()*V4_Eta_HeliEtaPipPim.Z()-V4_Pip_HeliEtaPipPim.Z()*V4_Eta_HeliEtaPipPim.Y(),
			V4_Pip_HeliEtaPipPim.Z()*V4_Eta_HeliEtaPipPim.X()-V4_Pip_HeliEtaPipPim.X()*V4_Eta_HeliEtaPipPim.Z(),
			V4_Pip_HeliEtaPipPim.X()*V4_Eta_HeliEtaPipPim.Y()-V4_Pip_HeliEtaPipPim.Y()*V4_Eta_HeliEtaPipPim.X());

		EvtData* evtData=new EvtData();
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Psi] = V4_psi;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPipPim_HeliPsi] = V4_EtaPipPim_HeliPsi;

		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi] = V4_PipPim_HeliPsi;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi] = V4_EtaPip_HeliPsi;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi] = V4_EtaPim_HeliPsi;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_gamma_HeliPsi] = V4_gamma_HeliPsi;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Eta_HeliPsi] = V4_Eta_HeliPsi;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pip_HeliPsi] = V4_Pip_HeliPsi;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pim_HeliPsi] = V4_Pim_HeliPsi;

		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Eta_HeliEtaPipPim] = V4_Eta_HeliEtaPipPim;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pip_HeliEtaPipPim] = V4_Pip_HeliEtaPipPim;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pim_HeliEtaPipPim] = V4_Pim_HeliEtaPipPim;

		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliEtaPipPim] = V4_EtaPip_HeliEtaPipPim;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliEtaPipPim] = V4_EtaPim_HeliEtaPipPim;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliEtaPipPim] = V4_PipPim_HeliEtaPipPim;

		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pip_HeliPipPim] = V4_Pip_HeliPipPim;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pim_HeliPipPim] = V4_Pim_HeliPipPim;

		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pip_HeliEtaPip] = V4_Pip_HeliEtaPip;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pim_HeliEtaPim] = V4_Pim_HeliEtaPim;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Eta_HeliEtaPip] = V4_Eta_HeliEtaPip;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_Eta_HeliEtaPim] = V4_Eta_HeliEtaPim;

		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_normPipPimDecHeliEtaPipPim] = V4_normPipPimDecHeliEtaPipPim;
		evtData->FourVecs[enumJpsiGamEtaPiPiData::V4_normEtaPipDecHeliEtaPipPim] = V4_normEtaPipDecHeliEtaPipPim;

		// calculate and store WignerD functions for Psi -> X gamma (JPC =0-=, 0++, 2++)
		Spin jPsi=1;
		for (Spin M=-1; M<=1; M=M+2){
			for (Spin lam=-1; lam<=1; lam++){
				evtData->WignerDs[enumJpsiGamEtaPiPiData::Df_Psi][jPsi][M][lam]=Wigner_D(0.,V4_EtaPipPim_HeliPsi.Theta(),0,jPsi,M,lam);
			}
		}

		/*		//WignerD function for 2+ -> phi phi
		Spin jTensor =2;
		for(Spin M=-jTensor; M<=jTensor; M++){
			for (Spin lam=-2; lam<=2; lam++){
				evtData->WignerDs[enumJpsiGamEtaPiPiData::Df_Spin2][jTensor][M][lam] = Wigner_D(V4_KsKl_HeliKsKlKpKm.Phi(),V4_KsKl_HeliKsKlKpKm.Theta(), 0,jTensor,M,lam);
			}
		}

		//WignerD function for 1+ -> phi phi
		Spin jAxialVec =1;
		for(Spin M=-jAxialVec; M<=jAxialVec; M++){
			for (Spin lam=-1; lam<=1; lam++){
				evtData->WignerDs[enumJpsiGamEtaPiPiData::Df_Spin2][jAxialVec][M][lam] = Wigner_D(V4_KsKl_HeliKsKlKpKm.Phi(),V4_KsKl_HeliKsKlKpKm.Theta(), 0,jAxialVec,M,lam); //use Df_Spin2 (!) to be able to use same function for fJ(J=1,2) decay amplitudes
			}
		}

		//WignerD function for 0-/0+ -> phi phi
		Spin jScalar =0;
		Spin M =0;
		Spin lam =0;
		evtData->WignerDs[enumJpsiGamEtaPiPiData::Df_Spin0][jScalar][M][lam] = Wigner_D(V4_KsKl_HeliKsKlKpKm.Phi(),V4_KsKl_HeliKsKlKpKm.Theta(), 0,jScalar,M,lam);


		//WignerD function for phi -> K+ K- and phi -> KS KL
		Spin phiSpin=1;
		for(Spin M=-phiSpin; M<=phiSpin; M++){
			Spin lam=0;
			evtData->WignerDs[enumJpsiGamEtaPiPiData::Df_KsKl][phiSpin][M][lam] = Wigner_D(V4_Ks_HeliKsKl.Phi(),V4_Ks_HeliKsKl.Theta(), 0,phiSpin,M,lam);
			evtData->WignerDs[enumJpsiGamEtaPiPiData::Df_KpKm][phiSpin][M][lam] = Wigner_D(V4_Kp_HeliKpKm.Phi(),V4_Kp_HeliKpKm.Theta(), 0,phiSpin,M,lam);
		}
		*/

		evtData->evtWeight=1.;
		theEvtList.push_back(evtData);

		++evtCount;
	}
}


