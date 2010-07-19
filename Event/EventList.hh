#ifndef EVENTLIST_HH
#define EVENTLIST_HH

#include <vector>

class Event;
class Particle;
class ParticleTable;

class EventList
{
public:
  EventList();
  ~EventList();

  void add(Event*);

  Event* nextEvent();
  void rewind();

  void setParticleType(int num, Particle* partRef);
  bool findParticleTypes(ParticleTable& pdtTable);

  Particle* particle(unsigned int num);

  int size();

private:
  std::vector<Event*> eventList;
  std::vector<Event*>::const_iterator currentEvent;
  std::vector<Particle*> particleRefs;
  
};

#endif
