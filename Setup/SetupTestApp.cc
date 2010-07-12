#include "Setup/PwaEnv.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "qft++/topincludes/tensor.hh"

#include "ErrLogger/ErrLogger.hh"

#include <iostream>
#include <string>
#include <cstdlib>

int main()
{
  //  ErrLineLog myLogger(ErrLog::debugging);
  
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR"); 
  std::string setupFile(theSourcePath + "/Setup/test.setup");
  PwaEnv::instance().setup(setupFile);

  ParticleTable* pTable = PwaEnv::instance().particleTable();
  if (0 == pTable) {
    Alert << "getting ParticleTable failed" << endmsg;
    exit(1);
  }
  pTable->print(std::cout);

  EventList* eventList = PwaEnv::instance().beamEventList();
  if (0 == eventList) {
    Alert << "getting beam EventList failed" << endmsg;
    exit(1);
  }

  Info << "Input file has " << eventList->size() << " events. Each event has "
       <<  eventList->nextEvent()->size() << " final state particles.\n" << endmsg;
  eventList->rewind();

  EventList* mcEventList = PwaEnv::instance().mcEventList();
  if (0 == mcEventList) {
    Alert << "getting MC EventList failed" << endmsg;
    exit(1);
  }

  Info << "MC Input file has " << mcEventList->size() << " events. Each event has "
       <<  mcEventList->nextEvent()->size() << " final state particles.\n" << endmsg;
  mcEventList->rewind();

  Event* anEvent;
  int evtCount = 0;
  Info << "======== beam events ========" << endmsg;
  while ((anEvent = eventList->nextEvent()) != 0 && evtCount < 20) {
    Info << "\n" 
	 << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
	 << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
	 << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
	 << endmsg;
    ++evtCount;
  }

  evtCount = 0;
  Info << "======== MC events ========" << endmsg;
  while ((anEvent = mcEventList->nextEvent()) != 0 && evtCount < 20) {
    Info << "\n" 
	 << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
	 << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
	 << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
	 << endmsg;
    ++evtCount;
  }


  return 0;
}
