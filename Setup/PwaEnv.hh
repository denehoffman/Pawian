#pragma once

#include <string>

class ParticleTable;
class EventList;
class DecayTree;

class PwaEnv
{
public:
  static PwaEnv& instance();
  ~PwaEnv();

  bool setup(std::string& setupFileName);

  ParticleTable* particleTable();
  EventList* beamEventList();
  EventList* mcEventList();

private:
  PwaEnv();

  ParticleTable* theParticleTable;
  DecayTree* theDecayTree;
  EventList* theBeamEventList;
  EventList* theMcEventList;
};
