#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh"

#include "ErrLogger/ErrLineLog.hh"

#include <iostream>
#include <string>
#include <cstdlib>

int main()
{

  ErrLineLog thisLogger(ErrLog::debugging);

  ParticleTable pTable;
  PdtParser parser;
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR"); 
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");

  if (parser.parse(pdtFile, pTable)) {
    pTable.print(std::cout);

    std::string name("omega");
    Particle* omega = pTable.particle(name);
    if (0 != omega)
      omega->print(std::cout);
    else
      ErrMsg(warning) << name << " not found" << endmsg;

    name = std::string("oohps");
    Particle* dummy = pTable.particle(name);
    if (0 != dummy)
      dummy->print(std::cout);
    else
      ErrMsg(warning) << name << " not found" << endmsg;

    // try clone error handling
    pTable.clone(std::string("omegaNew"), std::string("omeGa")); // syntax here is clone(new, old);
    pTable.clone(std::string("omegaNew"), std::string("omega"));
    pTable.particle(std::string("omegaNew"))->print(std::cout);

    pTable.clone(std::string("omegaNew"), std::string("omega"));

    // try to modify mass and width
    pTable.modifyMass(std::string("omegaNew"), 1.42);
    pTable.modifyWidth(std::string("omegaNew"), 0.13);
    pTable.particle(std::string("omegaNew"))->print(std::cout);

  } else {
    ErrMsg(error) << "Error: could not parse " << pdtFile << endmsg;
    exit(1);
  }

  return 0;
}
