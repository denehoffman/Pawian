#pragma once

#include <string>
#include "Particle/ParticleTable.hh"

struct ParticleData;

class PdtParser
{
public:

  PdtParser();
  ~PdtParser();

  bool parse(std::string& fileName, ParticleTable& table);
  bool parse(std::string::const_iterator begin,
	     std::string::const_iterator end,
	     ParticleData& pData);
};
