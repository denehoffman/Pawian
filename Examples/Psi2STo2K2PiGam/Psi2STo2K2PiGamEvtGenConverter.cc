#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtGenConverter.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Event/Event.hh"
#include "Event/EventList.hh"
#include "Event/HepMCEventList.hh"
#include "ErrLogger/ErrLogger.hh"
#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"

Psi2STo2K2PiGamEvtGenConverter* Psi2STo2K2PiGamEvtGenConverter::_instance=0;

Psi2STo2K2PiGamEvtGenConverter::Psi2STo2K2PiGamEvtGenConverter() :
  _psi2STo2K2PiGamEvtList(new Psi2STo2K2PiGamEvtList())
{
}

Psi2STo2K2PiGamEvtGenConverter::~Psi2STo2K2PiGamEvtGenConverter()
{
  delete _psi2STo2K2PiGamEvtList;
  delete _instance;
  _instance=0;
}


Psi2STo2K2PiGamEvtGenConverter* Psi2STo2K2PiGamEvtGenConverter::instance()
{
  if (_instance==0) _instance=new Psi2STo2K2PiGamEvtGenConverter();
  return _instance;
}

Event* Psi2STo2K2PiGamEvtGenConverter::convertHepMCtoEvent(HepMC::GenEvent* HepMCevt){

  Event* result=0;
  
  if(0==HepMCevt){
    Warning << "HepMCevt is 0 pointer; cannot convert HepMC event to PAWIAN Event; return 0 pointer!!!"
	    << endmsg;
    return result;
  }

  std::vector<HepMC::GenParticle*> KplusHepMC;
  std::vector<HepMC::GenParticle*> KminusHepMC;
  std::vector<HepMC::GenParticle*> PhotonRadHepMCVec;
  std::vector<HepMC::GenParticle*> Pi0HepMCVec;

  for ( HepMC::GenEvent::particle_iterator p = HepMCevt->particles_begin();
	p != HepMCevt->particles_end(); ++p ) {
    if ((*p)->pdg_id() == 22) PhotonRadHepMCVec.push_back(*p); //photons
    else if ((*p)->pdg_id() == 111) Pi0HepMCVec.push_back(*p); //pi0s
    else if ((*p)->pdg_id() == -321) KminusHepMC.push_back(*p); //K-s
    else if ((*p)->pdg_id() == 321) KplusHepMC.push_back(*p); //K-s
  }
  
  if ( PhotonRadHepMCVec.size()!=1 ||  Pi0HepMCVec.size()!=2 || KplusHepMC.size()!=1 || KminusHepMC.size()!=1){
    Alert << "extracted final state particles is not consistent with K+ K- 2pi0 gamma final state!!!" << endmsg;
    exit(0);  
  }

  int theEvtNo=HepMCevt->event_number();
  Event* newEvent = new Event(theEvtNo);  
  //now fill the event list: order: gamma, K+, K-, pi0, pi0
  //photon
  const HepMC::FourVector radGam4Vec=PhotonRadHepMCVec[0]->momentum();
  newEvent->addParticle(radGam4Vec.e(), radGam4Vec.px(), radGam4Vec.py(), radGam4Vec.pz());
  
  //K+   
  const HepMC::FourVector Kplus4Vec=KplusHepMC[0]->momentum();
  newEvent->addParticle(Kplus4Vec.e(), Kplus4Vec.px(),Kplus4Vec.py(), Kplus4Vec.pz());
  
  //K-   
  const HepMC::FourVector Kminus4Vec=KminusHepMC[0]->momentum();
  newEvent->addParticle(Kminus4Vec.e(), Kminus4Vec.px(),Kminus4Vec.py(), Kminus4Vec.pz());
  
  //pi0's 
  for(int i=0; i<2; i++){
    const HepMC::FourVector pi4Vec=Pi0HepMCVec[i]->momentum();
    newEvent->addParticle(pi4Vec.e(), pi4Vec.px(), pi4Vec.py(), pi4Vec.pz());
  }


  return newEvent;
}

Psi2STo2K2PiGamEvtData* Psi2STo2K2PiGamEvtGenConverter::convertHepMCtoDataEvt(HepMC::GenEvent* genEvent){

  Event* tmpEvent=convertHepMCtoEvent(genEvent);
  Psi2STo2K2PiGamEvtData* theEvtData = _psi2STo2K2PiGamEvtList->fillEvtData(tmpEvent);
  delete tmpEvent;
  return theEvtData;
}

bool Psi2STo2K2PiGamEvtGenConverter::fillEventList(boost::shared_ptr<HepMCEventList> hepMCEventListPtr, EventList& evtList){

  hepMCEventListPtr->rewind();
  HepMC::GenEvent* genEvent;
  int evtCount = 0;
  while ( (genEvent = hepMCEventListPtr->nextEvent()) !=0 ){
    genEvent->print(std::cout);
    Event* newEvent=convertHepMCtoEvent(genEvent);
    
    evtList.add(newEvent);
    
    evtCount++;
  }

  hepMCEventListPtr->rewind();  
  Info << "HepMCEventList contains " << evtCount << " events" << endmsg;
  return true; // success
}

bool Psi2STo2K2PiGamEvtGenConverter::fillEventMap(boost::shared_ptr<HepMCEventList> hepMCEventListPtr, std::map<Event*, HepMC::GenEvent*, pawian::Collection::PtrLess>& theEvtMap){

  hepMCEventListPtr->rewind();
  HepMC::GenEvent* genEvent;
  int evtCount = 0;
  while ( (genEvent = hepMCEventListPtr->nextEvent()) !=0 ){
    genEvent->print(std::cout);
    Event* newEvent=convertHepMCtoEvent(genEvent);
    theEvtMap[newEvent]=genEvent;
    
    evtCount++;
  }

  hepMCEventListPtr->rewind();  
  return true; // success

}

bool Psi2STo2K2PiGamEvtGenConverter::fillEventMap(boost::shared_ptr<HepMCEventList> hepMCEventListPtr, std::map<Psi2STo2K2PiGamEvtData*, HepMC::GenEvent*, pawian::Collection::PtrLess>& theEvtMap){

  hepMCEventListPtr->rewind();
  HepMC::GenEvent* genEvent;
  int evtCount=0;
  while ( (genEvent = hepMCEventListPtr->nextEvent()) !=0 ){
    genEvent->print(std::cout);
    Event* newEvent=convertHepMCtoEvent(genEvent);
    Psi2STo2K2PiGamEvtData* theEvtData=_psi2STo2K2PiGamEvtList->fillEvtData(newEvent, evtCount); 
    theEvtMap[theEvtData]=genEvent;
    evtCount++;
    delete newEvent;
  }
  return true;
}
