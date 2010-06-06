#include "Setup/PwaEnv.hh"

#include "Particle/ParticleTable.hh"
#include "Particle/PdtParser.hh"
#include "Setup/SetupParser.hh"
#include "Event/EventList.hh"
#include "ErrLogger/ErrLineLog.hh"
#include "Event/CBElsaReader.hh"

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
    ErrMsg(fatal) << "can not parse particle table " << pdtFile << endmsg;
  }

  SetupParser setupParser;
  if (!setupParser.parse(setupFileName, theParticleTable)) {
    ErrMsg(fatal) << "can not parse setup " << setupFileName << endmsg;
  }

  if (setupParser.setup()->beamInput.size() == 0)
    ErrMsg(fatal) << "no beam data defined in setup" << endmsg;

  CBElsaReader eventReader(setupParser.setup()->beamInput, 3, 1);
  theBeamEventList = new EventList();
  eventReader.fillAll(*theBeamEventList);

  if (setupParser.setup()->mcInput.size() == 0)
    ErrMsg(fatal) << "no MC data defined in setup" << endmsg;

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
