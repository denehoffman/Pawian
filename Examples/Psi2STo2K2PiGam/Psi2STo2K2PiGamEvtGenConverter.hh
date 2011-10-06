#ifndef Psi2STo2K2PiGamEvtGenConverter_HH
#define Psi2STo2K2PiGamEvtGenConverter_HH

#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "HepMC/GenEvent.h"
#include "Utils/PawianCollectionUtils.hh"
#include <map>
#include <boost/shared_ptr.hpp>

using Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData;

class Event;
class EventList;
class HepMCEventList;
class Psi2STo2K2PiGamEvtList;

class Psi2STo2K2PiGamEvtGenConverter
{

public:
  static Psi2STo2K2PiGamEvtGenConverter* instance();
  ~Psi2STo2K2PiGamEvtGenConverter();

  Event* convertHepMCtoEvent(HepMC::GenEvent*);
  Psi2STo2K2PiGamEvtData* convertHepMCtoDataEvt(HepMC::GenEvent*);
  bool fillEventList(boost::shared_ptr<HepMCEventList> hepMCEventListPtr, EventList& evtList);
  bool fillEventMap(boost::shared_ptr<HepMCEventList> hepMCEventListPtr, std::map<Event*, HepMC::GenEvent*, pawian::Collection::PtrLess>& theEvtMap);
  bool fillEventMap(boost::shared_ptr<HepMCEventList> hepMCEventListPtr, std::map<Psi2STo2K2PiGamEvtData*, HepMC::GenEvent*, pawian::Collection::PtrLess>& theEvtMap);

private:
  Psi2STo2K2PiGamEvtGenConverter();

  static Psi2STo2K2PiGamEvtGenConverter* _instance;
  Psi2STo2K2PiGamEvtList* _psi2STo2K2PiGamEvtList; 
};

#endif
