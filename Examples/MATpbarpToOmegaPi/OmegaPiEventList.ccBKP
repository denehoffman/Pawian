#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/MATpbarpToOmegaPi/OmegaPiEventList.hh"
#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"

#include "Examples/MATpbarpToOmegaPi/serSpin.hh"
#include "Examples/MATpbarpToOmegaPi/serVec4.hh"

OmegaPiEventList::OmegaPiEventList(EventList& evtListData, EventList& evtListMc, unsigned jmax, unsigned pbarmom):
  _jmax(jmax),
  _pbarmom(pbarmom)
{
  if (_jmax<0) { Alert << "_jmax < 0 is not allowed!!!" << endmsg; exit(1); }
  read4Vecs(evtListData, _dataList);
  read4Vecs(evtListMc, _mcList);
}

OmegaPiEventList::OmegaPiEventList()
{
}

OmegaPiEventList::~OmegaPiEventList()
{
}

void OmegaPiEventList::read4Vecs(EventList& evtList, std::vector<OmPiEvtData*>& omPiEvtList)
{

  Event* anEvent;
  int evtCount = 0;
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

    serVec4<float> omega_cm_4V(omega_4V);
    omega_cm_4V.Boost(cm_4V);

    serVec4<float> piRec_cm_4V(piRec_4V);
    piRec_cm_4V.Boost(cm_4V);

    serVec4<float> pi0FromOmega4V=*(anEvent->p4(2));
    if ( fabs(pi0FromOmega4V.M()-0.13497)>0.01 ) {
      Alert <<"the third particle is not the pi0 from the omega decay" << endmsg;
      exit(1);
    }

    serVec4<float> pi0HeliOmega4V = helicityVec(cm_4V, omega_4V, pi0FromOmega4V);

    OmPiEvtData* theOmPiEvtData=new OmPiEvtData();
   // for(int i=0; i<4; i++){
   //   theOmPiEvtData->cm_4Vec.push_back(cm_4V[i]);  //theOmPiEvtData->cm_4Vec=cm_4V;
   //   theOmPiEvtData->omegaHeliCm4Vec.push_back(omega_cm_4V[i]);  //theOmPiEvtData->omegaHeliCm4Vec=omega_cm_4V;
   //   theOmPiEvtData->pi0RecHeliCm4Vec.push_back(piRec_cm_4V[i]);  //theOmPiEvtData->pi0RecHeliCm4Vec=piRec_cm_4V;
   //   theOmPiEvtData->pi0HeliOmega4Vec.push_back(pi0HeliOmega4V[i]);  //theOmPiEvtData->pi0HeliOmega4Vec=pi0HeliOmega4V;
   // }

    for (Spin j=0; j<=_jmax; j++){
      for (Spin M=-1; M<=1; M++){
	if (fabs(M)>j) continue;
	for (Spin lam=-1; lam<=1; lam++){
 	  if (fabs(lam)>j) continue;
	  theOmPiEvtData->Dfp[serSpin(j)][serSpin(M)][serSpin(lam)]=Wigner_D(omega_cm_4V.Phi(),omega_cm_4V.Theta(),0,j,M,lam);
	}
      }
    }

    Spin omegaSpin=1;
    myWigner_D(omegaSpin, pi0HeliOmega4V.Phi(), pi0HeliOmega4V.Theta(),0, theOmPiEvtData->Dfd);
    
    omPiEvtList.push_back(theOmPiEvtData);

    ++evtCount;
  }
 
}

void OmegaPiEventList::myWigner_D(const Spin &__jmax,double __alpha,double __beta,double __gamma,
	      map<serSpin,map<serSpin,map<serSpin,complex<double> > > > &__D){

  complex<double> i(0.,1.);
  map<Spin,map<Spin,map<Spin,double> > > d;
  Wigner_d(__jmax,__beta,d);
  Spin jmin;
  if(d.find(0) != d.end()) jmin = 0;
  else jmin = 1/2.;
  
  for(Spin j = jmin; j <= __jmax; j++){
    for(Spin m = -j; m <= j; m++){
      for(Spin n = -j; n <= j; n++) 
	__D[serSpin(j)][serSpin(m)][serSpin(n)] = exp(-i*(m*__alpha + n*__gamma))*d[j][m][n];
    }
  }
}

