#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


OmegaPiEventList::OmegaPiEventList(EventList& evtListData, EventList& evtListMc, unsigned jmax, unsigned pbarmom):
  _jmax(jmax),
  _pbarmom(pbarmom)
{
  if (_jmax<0) { Alert << "_jmax < 0 is not allowed!!!" << endmsg; exit(1); }
  read4Vecs(evtListData, _dataList);
  read4Vecs(evtListMc, _mcList);
}



OmegaPiEventList::~OmegaPiEventList()
{
}

void OmegaPiEventList::read4Vecs(EventList& evtList, std::vector<OmPiEvtData*>& omPiEvtList)
{

  Event* anEvent;
  int evtCount = 0;
  bool bPbarmomOnce=false;
  while ((anEvent = evtList.nextEvent()) != 0){
    if ( evtCount%1000 == 0 ) Info << "4vec calculation for event " << evtCount << endmsg;
    Vector4<float>  cm_4V(*(anEvent->p4(0))+*(anEvent->p4(1)));

    Vector4<float> omega_4V;
    Vector4<float> piRec_4V;


    for(int i=0; i<2; ++i){
      if ( fabs(anEvent->p4(i)->Mass()-0.78195)<0.01 ) omega_4V=*(anEvent->p4(i));
      else if ( fabs(anEvent->p4(i)->Mass()-0.13497)<0.01 ) piRec_4V=*(anEvent->p4(i));
      else {
	Alert <<"this is neither an omega nor a pi0 particle!!!" << endmsg;
	exit(1);
      }
    }

    Vector4<float> omega_cm_4V(omega_4V);
    omega_cm_4V.Boost(cm_4V);

    Vector4<float> piRec_cm_4V(piRec_4V);
    piRec_cm_4V.Boost(cm_4V);
    


  if (!bPbarmomOnce)
   {
     float E=sqrt(pow(938.272,2)+pow(double(_pbarmom),2));
     Vector4<float> pbarmom;
     pbarmom.SetP4(E,0.,0.,_pbarmom);
     pbarmom.Boost(cm_4V);
     cout << endl << "pbarmom-cmc =" << pbarmom.P() << " MeV/c" << endl;
     bPbarmomOnce=true;
   }
    
    Vector4<float> pi0FromOmega4V=*(anEvent->p4(2));
    if ( fabs(pi0FromOmega4V.M()-0.13497)>0.01 ) {
      Alert <<"the third particle is not the pi0 from the omega decay" << endmsg;
      exit(1);
    }

    Vector4<float> pi0HeliOmega4V = helicityVec(cm_4V, omega_4V, pi0FromOmega4V);

    Vector4<float> piDec_cm_4V(pi0FromOmega4V);
    piDec_cm_4V.Boost(cm_4V);


    OmPiEvtData* theOmPiEvtData=new OmPiEvtData();
    theOmPiEvtData->cm_4Vec=cm_4V;
    theOmPiEvtData->omegaHeliCm4Vec=omega_cm_4V;
    theOmPiEvtData->pi0RecHeliCm4Vec=piRec_cm_4V;
    theOmPiEvtData->pi0DecHeliCm4Vec=piDec_cm_4V;
    theOmPiEvtData->pi0HeliOmega4Vec=pi0HeliOmega4V;
    theOmPiEvtData->cosPi0HeliOmega4Vec=costDecHeli(cm_4V, omega_4V, pi0FromOmega4V);

    for (Spin j=0; j<=_jmax; j++){
      for (Spin M=-1; M<=1; M++){
	if (fabs(M)>j) continue;
	for (Spin lam=-1; lam<=1; lam++){
 	  if (fabs(lam)>j) continue;
	  theOmPiEvtData->Dfp[j][M][lam]=Wigner_D(0.,omega_cm_4V.Theta(),0,j,M,lam);
	}
      }
    }

    Spin omegaSpin=1;
    Wigner_D(omegaSpin, pi0HeliOmega4V.Phi(), pi0HeliOmega4V.Theta(),0, theOmPiEvtData->Dfd);
    
    omPiEvtList.push_back(theOmPiEvtData);

    ++evtCount;
  }
 
}


