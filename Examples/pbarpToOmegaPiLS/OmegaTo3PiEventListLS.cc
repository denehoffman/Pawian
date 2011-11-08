#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPiLS/OmegaTo3PiEventListLS.hh"
#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


OmegaTo3PiEventListLS::OmegaTo3PiEventListLS(EventList& evtListData, EventList& evtListMc, unsigned jmax, unsigned pbarmom):
  AbsOmegaPiEventListLS(evtListData, evtListMc, jmax, pbarmom)
{
  read4Vecs(evtListData, _dataList);
  read4Vecs(evtListMc, _mcList);
}




OmegaTo3PiEventListLS::~OmegaTo3PiEventListLS()
{
}

void OmegaTo3PiEventListLS::read4Vecs(EventList& evtList, std::vector<OmPiEvtDataLS*>& omPiEvtList)
{

  Event* anEvent;
  int evtCount = 0;
  bool bPbarmomOnce=false;
  while ((anEvent = evtList.nextEvent()) != 0){
    if ( evtCount%1000 == 0 ) Info << "4vec calculation for event " << evtCount << endmsg;
    Vector4<float>  cm_4V( *(anEvent->p4(0))+ *(anEvent->p4(1)) + *(anEvent->p4(2)) + *(anEvent->p4(3)));

    Vector4<float> pi0Omega_4V;
    Vector4<float> pipOmega_4V;
    Vector4<float> pimOmega_4V;
    Vector4<float> pi0Rec_4V;

    pi0Rec_4V=*(anEvent->p4(0));
    pi0Omega_4V=*(anEvent->p4(1));
    pipOmega_4V=*(anEvent->p4(2));
    pimOmega_4V=*(anEvent->p4(3));

 
    Vector4<float> omega_4V(pi0Omega_4V+pipOmega_4V+pimOmega_4V);

    Vector4<float> omega_cm_4V(omega_4V);
    omega_cm_4V.Boost(cm_4V);

    Vector4<float> pi0Rec_cm_4V(pi0Rec_4V);
    pi0Rec_cm_4V.Boost(cm_4V);
    


  if (!bPbarmomOnce)
   {
     float E=sqrt(pow(938.272,2)+pow(double(_pbarmom),2));
     Vector4<float> pbarmom;
     pbarmom.SetP4(E,0.,0.,_pbarmom);
     pbarmom.Boost(cm_4V);
     cout << endl << "pbarmom-cmc =" << pbarmom.P() << " MeV/c" << endl;
     bPbarmomOnce=true;
   }
    

    Vector4<float> omega_4V2(pi0Rec_4V+pipOmega_4V+pimOmega_4V);
    Vector4<float> omega_cm_4V2(omega_4V2);
    omega_cm_4V2.Boost(cm_4V);

    Vector4<float> pi0HeliOmega4V = helicityVec(cm_4V, omega_4V, pi0Omega_4V);
    Vector4<float> pipHeliOmega4V = helicityVec(cm_4V, omega_4V, pipOmega_4V);
    Vector4<float> pimHeliOmega4V = helicityVec(cm_4V, omega_4V, pimOmega_4V);

    Vector4<float> pi0HeliOmega4V2 = helicityVec(cm_4V, omega_4V2, pi0Rec_4V);
    Vector4<float> pipHeliOmega4V2 = helicityVec(cm_4V, omega_4V2, pipOmega_4V);
    Vector4<float> pimHeliOmega4V2 = helicityVec(cm_4V, omega_4V2, pimOmega_4V);

    Vector4<float> normOmegaDecHeliOmega_4V(0.5*(pi0HeliOmega4V.T()+pipHeliOmega4V.T()+pimHeliOmega4V.T()),
					    pimHeliOmega4V.Y()*pipHeliOmega4V.Z()-pimHeliOmega4V.Z()*pipHeliOmega4V.Y(),
					    pimHeliOmega4V.Z()*pipHeliOmega4V.X()-pimHeliOmega4V.X()*pipHeliOmega4V.Z(),		
					    pimHeliOmega4V.X()*pipHeliOmega4V.Y()-pimHeliOmega4V.Y()*pipHeliOmega4V.X());			

    Vector4<float> normOmegaDecHeliOmega_4V2(0.5*(pi0HeliOmega4V2.T()+pipHeliOmega4V2.T()+pimHeliOmega4V2.T()),
					    pimHeliOmega4V2.Y()*pipHeliOmega4V2.Z()-pimHeliOmega4V2.Z()*pipHeliOmega4V2.Y(),
					    pimHeliOmega4V2.Z()*pipHeliOmega4V2.X()-pimHeliOmega4V2.X()*pipHeliOmega4V2.Z(),		
					    pimHeliOmega4V2.X()*pipHeliOmega4V2.Y()-pimHeliOmega4V2.Y()*pipHeliOmega4V2.X());		    


    float theQ=pipHeliOmega4V.E()-pipHeliOmega4V.M()+pimHeliOmega4V.E()-pimHeliOmega4V.M()+pi0HeliOmega4V.E()-pi0HeliOmega4V.M();
    float lambdaNorm=theQ*theQ*(theQ*theQ/108.+pimHeliOmega4V.M()*theQ/9.+pimHeliOmega4V.M()*pimHeliOmega4V.M()/3.);
    float lambda=normOmegaDecHeliOmega_4V.P()*normOmegaDecHeliOmega_4V.P()/lambdaNorm;

    OmPiEvtDataLS* theOmPiEvtData=new OmPiEvtDataLS();
    theOmPiEvtData->cm_4Vec=cm_4V;
    theOmPiEvtData->omegaHeliCm4Vec=omega_cm_4V;
    theOmPiEvtData->omegaHeliCm4Vec2=omega_cm_4V2;

    theOmPiEvtData->pi0RecHeliCm4Vec=pi0Rec_cm_4V;


    theOmPiEvtData->pi0HeliOmega4Vec=normOmegaDecHeliOmega_4V;
    theOmPiEvtData->pi0HeliOmega4Vec2=normOmegaDecHeliOmega_4V2;

    theOmPiEvtData->cosPi0HeliOmega4Vec=normOmegaDecHeliOmega_4V.CosTheta();

    theOmPiEvtData->pipHeliOmega4Vec=pipHeliOmega4V;
    theOmPiEvtData->pimHeliOmega4Vec=pimHeliOmega4V;
    theOmPiEvtData->pippimCrossHeliOmega4Vec=normOmegaDecHeliOmega_4V;
    theOmPiEvtData->lambda=lambda;
    theOmPiEvtData->eventWeight = anEvent->Weight();

    for (Spin j=0; j<=_jmax; j++){
      for (Spin M=-1; M<=1; M++){
	if (fabs(M)>j) continue;
        for (Spin lam=-1; lam<=1; lam++){
          if (fabs(lam)>j) continue;
	  theOmPiEvtData->Dfp[j][M][lam]=Wigner_D(0.,omega_cm_4V.Theta(),0,j,M,lam);
	  theOmPiEvtData->Dfp2[j][M][lam]=Wigner_D(0.,omega_cm_4V2.Theta(),0,j,M,lam);
	}
      }
    }

    Spin omegaSpin=1;
    Wigner_D(omegaSpin, normOmegaDecHeliOmega_4V.Phi(), normOmegaDecHeliOmega_4V.Theta(),0, theOmPiEvtData->Dfd);
    Wigner_D(omegaSpin, normOmegaDecHeliOmega_4V2.Phi(), normOmegaDecHeliOmega_4V2.Theta(),0, theOmPiEvtData->Dfd2);
    
  
   
    omPiEvtList.push_back(theOmPiEvtData);

    ++evtCount;
  }
 
}


