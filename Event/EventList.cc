#include "Event/EventList.hh"

#include "Event/Event.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "ErrLogger/ErrLineLog.hh"

EventList::EventList()
{
}

EventList::~EventList()
{
  for (currentEvent = eventList.begin();
       currentEvent != eventList.end();
       ++currentEvent)
    delete *currentEvent;
}

void EventList::add(Event* newEvent)
{
  if (0 != newEvent)
    eventList.push_back(newEvent);
  else
    ErrMsg(fatal) << "can not add 0 pointer to event list" << endmsg;

  return;
}

Event* EventList::nextEvent()
{
  if (currentEvent != eventList.end()) {
    Event* result = *currentEvent;
    ++currentEvent;
    return result;
  } else
    return 0;
}
 
void EventList::rewind()
{
  currentEvent = eventList.begin();
  return;
}

void EventList::setParticleType(int num, Particle* partRef)
{
  particleRefs[num] = partRef;
  return;
}

Particle* EventList::particle(int num)
{
  if (num >= 0 && num < particleRefs.size())
    return particleRefs[num];

  return 0;
}

int EventList::size()
{
  return eventList.size();
}

bool EventList::findParticleTypes(ParticleTable& pdtTable)
{
  if (eventList.size() == 0)
    ErrMsg(fatal) << "eventList is empty" << endmsg;

  bool result = true;
  int partsPerEvent = eventList[0]->size();
  int i;
  for (i=0; i<partsPerEvent; i++) {
    Particle* thisParticle;
    thisParticle = pdtTable.particle(eventList[0]->p4(i)->Mass());
    particleRefs.push_back(thisParticle);
    if (0 == thisParticle) {
      ErrMsg(warning) << "did not find a particle with mass " 
		      << eventList[0]->p4(i)->Mass() << endmsg;
      result = false;
    }
  }

  return result;
}
