#include "ErrLogger/ErrLogger.hh"
#include "Event/HepMCEventList.hh"
#include "HepMC/GenEvent.h"
#include <vector>
#include <string>

int main()
{
  ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
  std::vector<std::string> fileNames;

  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR"); 

  std::string hepMCFile(theSourcePath+"/Event/HepMCEvt.out");

  HepMCEventList hepMCEvtList(hepMCFile);
  hepMCEvtList.rewind();

  HepMC::GenEvent* genEvent;
  int evtCount = 0;
  while ( (genEvent = hepMCEvtList.nextEvent()) !=0 ){
    genEvent->print(std::cout);
    evtCount++;
  }

  Info << "HepMCEventList contains " << evtCount << " events" << endmsg;


  exit(0);
}
