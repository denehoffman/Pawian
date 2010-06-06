#ifndef PWAENV_HH
#define PWAENV_HH

#include <string>

class ParticleTable;
class EventList;

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
  EventList* theBeamEventList;
  EventList* theMcEventList;

};

#endif
