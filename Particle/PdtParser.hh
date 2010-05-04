#ifndef PDTPARSER_HH
#define PDTPARSER_HH

#include <string>
#include "ParticleTable.hh"


class PdtParser
{
public:

  PdtParser();
  ~PdtParser();

  bool parse(std::string& fileName, ParticleTable& table);

};

#endif
