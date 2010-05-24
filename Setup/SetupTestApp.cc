#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh"
#include "Setup/SetupParser.hh"

#include "ErrLogger/ErrLineLog.hh"

#include <iostream>
#include <string>
#include <cstdlib>

int main()
{
  ErrLineLog myLogger(ErrLog::debugging);

  ParticleTable pTable;
  PdtParser parser;
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR"); 
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");

  if (parser.parse(pdtFile, pTable)) {
    pTable.print(std::cout);

  } else {
    ErrMsg(fatal) << "Error: could not parse " << pdtFile << endmsg;
    exit(1);
  }

  SetupParser setupP;
  std::string setupFile(theSourcePath + "/Setup/test.setup");

  if (setupP.parse(setupFile)) {
    ErrMsg(trace) << "setup parsed successfully" << endmsg;
  } else {
    ErrMsg(fatal) << "Error: could not parse " << setupFile << endmsg;
    exit(1);
  }

  return 0;
}
