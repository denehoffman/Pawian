#include <iostream>
#include "Event/HepMCEventList.hh"
#include "ErrLogger/ErrLogger.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

HepMCEventList::HepMCEventList()
{
}

HepMCEventList::HepMCEventList(const std::vector<std::string>& files)
{
  fillFromFiles(files);
}

HepMCEventList::HepMCEventList(const std::string& file){
  fillFromFile(file);
}

HepMCEventList::~HepMCEventList()
{
  for (currentEvent = eventList.begin();
       currentEvent != eventList.end();
       ++currentEvent)
    delete *currentEvent;
}

void HepMCEventList::add(HepMC::GenEvent* newEvent)
{
  if (0 != newEvent)
    eventList.push_back(newEvent);
  else {
    Alert << "can not add 0 pointer to event list" << endmsg;
    exit(1);
  }

  return;
}

void HepMCEventList::removeEvents(unsigned int nBegin, unsigned int nEnd)
{
  if ( nBegin < nEnd || eventList.size()>nEnd) {
    Alert << "can not remove event no " << nBegin << " - " << nEnd 
	  << " from list" << endmsg; 
    exit(1);  
  }

  eventList.erase(eventList.begin()+nBegin, eventList.begin()+nEnd);
}

bool HepMCEventList::fillFromFile(const std::string& file){
  bool result=true;
  std::ifstream* currentStream= new std::ifstream(file.c_str(), std::ios::in);
  if (!currentStream) {
    Alert << "can not open " << file << endmsg;
    result=false;
    exit(1);
  }

  HepMC::IO_GenEvent* hepIOin=new HepMC::IO_GenEvent(*currentStream);
  
  hepIOin->use_input_units( HepMC::Units::GEV, HepMC::Units::MM );
  
  HepMC::GenEvent* evt = hepIOin->read_next_event();

  while ( evt ){
    //        if( evt->is_valid() ) {
//     evt->print(std::cout);
    //       }
    add(evt);
    (*hepIOin) >> evt;
  }

  delete hepIOin;
  delete currentStream;

  return result;
}

bool HepMCEventList::fillFromFiles(const std::vector<std::string>& files){

  bool result=true;
  std::vector<std::string>::const_iterator iter = files.begin();
  for (; iter != files.end(); ++iter){
    if (fillFromFile(*iter)==false) result=false;
  }
  return result;
}

HepMC::GenEvent* HepMCEventList::nextEvent()
{
  if (currentEvent != eventList.end()) {
    HepMC::GenEvent* result = *currentEvent;
    ++currentEvent;
    return result;
  } else
    return 0;
}
 
void HepMCEventList::rewind()
{
  currentEvent = eventList.begin();
  return;
}

int HepMCEventList::size()
{
  return eventList.size();
}

