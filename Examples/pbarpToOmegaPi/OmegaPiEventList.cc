#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLineLog.hh"


OmegaPiEventList::OmegaPiEventList(EventList& evtListData, EventList& evtListMc)
{
  read4Vecs(evtListData, _dataList);
  read4Vecs(evtListMc, _mcList);
}



OmegaPiEventList::~OmegaPiEventList()
{
}

void OmegaPiEventList::read4Vecs(EventList& evtList, std::vector<OmPiEvtData>& omPiEvtList)
{

  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = evtList.nextEvent()) != 0){
    if ( evtCount%1000 == 0 ) ErrMsg(routine) << "4vec calculation for event " << evtCount << endmsg;
    Vector4<float>  cm_4V(*(anEvent->p4(0))+*(anEvent->p4(1)));

    Vector4<float> omega_4V;
    Vector4<float> piRec_4V;

    for(int i=0; i<2; ++i){
      if ( fabs(anEvent->p4(i)->Mass()-0.78195)<0.01 ) omega_4V=*(anEvent->p4(i));
      else if ( fabs(anEvent->p4(i)->Mass()-0.13497)<0.01 ) piRec_4V=*(anEvent->p4(i));
      else {
	ErrMsg(fatal) <<"this is neither an omega nor a pi0 particle!!!" << endmsg;
      }
    }

    Vector4<float> omega_cm_4V(omega_4V);
    omega_cm_4V.Boost(cm_4V);

    Vector4<float> piRec_cm_4V(piRec_4V);
    piRec_cm_4V.Boost(cm_4V);

    Vector4<float> pi0FromOmega4V=*(anEvent->p4(2));
    if ( fabs(pi0FromOmega4V.M()-0.13497)>0.01 ) 
      ErrMsg(fatal) <<"the third particle is not the pi0 from the omega decay" << endmsg;

    Vector4<float> pi0HeliOmega4 = helicityVec(cm_4V, omega_4V, pi0FromOmega4V);

    OmPiEvtData theOmPiEvtData;
    theOmPiEvtData.cm_4Vec=cm_4V;
    theOmPiEvtData.omegaHeliCm4Vec=omega_cm_4V;
    theOmPiEvtData.pi0RecHeliCm4Vec=piRec_cm_4V;
    theOmPiEvtData.pi0HeliOmega4Vec=pi0HeliOmega4;

    omPiEvtList.push_back(theOmPiEvtData);

    ++evtCount;
  }
 
}


