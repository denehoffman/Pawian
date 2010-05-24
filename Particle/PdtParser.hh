#ifndef PDTPARSER_HH
#define PDTPARSER_HH

#include <string>
#include "Particle/ParticleTable.hh"


class PdtParser
{
public:

  PdtParser();
  ~PdtParser();

  bool parse(std::string& fileName, ParticleTable& table);

};

#endif
