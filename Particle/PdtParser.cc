#include "Particle/PdtParser.hh"

#include "Particle/ParticleGrammar.hh"

PdtParser::PdtParser()
{
}

PdtParser::~PdtParser()
{
}

bool PdtParser::parse(std::string& fileName, ParticleTable& table)
{
  ParticleData* pData;

  std::ifstream file;
  file.open(fileName.c_str());

  if (file.is_open()) {
    std::string str;
    while (!file.eof()) {
      std::getline(file, str);
      if (!file.eof()) {
	pData = new ParticleData;
	if (particleGrammar::parse_particle(str.begin(), str.end(), *pData)) { // success
	  Particle* newParticle = new Particle(*pData);
	  table.addParticle(newParticle);
	}
      }
    }
    file.close();
    return true; // success
  }
  else
    return false;
}

bool PdtParser::parse(std::string::const_iterator begin,
		      std::string::const_iterator end,
		      ParticleData& pData)
{
  return particleGrammar::parse_particle(begin, end, pData);
}
