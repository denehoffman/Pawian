#include "Setup/PwaEnv.hh"

#include "Particle/ParticleTable.hh"
#include "Particle/PdtParser.hh"
#include "Setup/SetupParser.hh"
#include "Event/EventList.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Event/CBElsaReader.hh"
#include "DecayTree/DecayTree.hh"

PwaEnv& PwaEnv::instance()
{
  static PwaEnv theEnv;
  return theEnv;
}

bool PwaEnv::setup(std::string& setupFileName)
{
  PdtParser pdtParser;
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR"); 
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");
  theParticleTable = new ParticleTable;

  if (!pdtParser.parse(pdtFile, *theParticleTable)) {
    Alert << "can not parse particle table " << pdtFile << endmsg;
    exit(1);
  }

  SetupParser setupParser;
  if (!setupParser.parse(setupFileName, theParticleTable)) {
    Alert << "can not parse setup " << setupFileName << endmsg;
    exit(1);
  }

  theDecayTree = new DecayTree(setupParser.edgeList());
  if (0 == theDecayTree) {
    Alert << "can not build decay tree" << endmsg;
    exit(1);
  }

  if (setupParser.setup()->beamInput.size() == 0) {
    Alert << "no beam data defined in setup" << endmsg;
    exit(1);
  }

  CBElsaReader eventReader(setupParser.setup()->beamInput, 3, 1);
  theBeamEventList = new EventList();
  eventReader.fillAll(*theBeamEventList);

  if (setupParser.setup()->mcInput.size() == 0) {
    Alert << "no MC data defined in setup" << endmsg;
    exit(1);
  }

  CBElsaReader mcEventReader(setupParser.setup()->mcInput, 3, 1);
  theMcEventList = new EventList();
  mcEventReader.fillAll(*theMcEventList);
  
}

ParticleTable* PwaEnv::particleTable()
{
  return theParticleTable;
}

EventList* PwaEnv::beamEventList()
{
  return theBeamEventList;
}

EventList* PwaEnv::mcEventList()
{
  return theMcEventList;
}

PwaEnv::PwaEnv()
{
}

PwaEnv::~PwaEnv()
{
  delete theParticleTable;
  delete theBeamEventList;
  delete theMcEventList;
}
