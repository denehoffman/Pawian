#include "Event/CBElsaReader.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/PdtParser.hh"
#include "ErrLogger/ErrLogger.hh"

#include "qft++/topincludes/tensor.hh"
#include <vector>
#include <string>

int main()
{
  ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
  std::vector<std::string> fileNames;

  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR"); 

  ParticleTable pTable;
  PdtParser parser;
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");
  if (!parser.parse(pdtFile, pTable)) {
    Alert << "Error: could not parse " << pdtFile << endmsg;
    exit(1);
  }

  std::string pipiomegaFile(theSourcePath+"/Event/710_1350.dat");
  fileNames.push_back(pipiomegaFile);

  CBElsaReader eventReader(fileNames, 3, 1);
  EventList pipiomegaEvents;
  eventReader.fillAll(pipiomegaEvents);

  if (!pipiomegaEvents.findParticleTypes(pTable))
    Warning << "could not find all particles" << endmsg;

  Info << "\nFile has " << pipiomegaEvents.size() << " events. Each event has "
       <<  pipiomegaEvents.nextEvent()->size() << " final state particles.\n" << endmsg;
  pipiomegaEvents.rewind();

  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = pipiomegaEvents.nextEvent()) != 0 && evtCount < 20) {
    Info << "\n" 
	 << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
	 << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
	 << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
	 << endmsg;
    ++evtCount;
  }

  exit(0);
}
