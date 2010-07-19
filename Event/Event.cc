#include "Event/Event.hh"
#include <vector>

Event::Event()
{}

Event::~Event()
{
  std::vector<EvtPartData*>::iterator iter;
  for (iter = particles.begin(); iter != particles.end(); ++iter)
    delete *iter;
}
  
Vector4<float>* Event::p4(unsigned int i)
{
  if (particles.size() > i)
    return particles[i]->vector4;
  else {
    Alert << "accessing " << i << "th of " 
	  << particles.size() << " particles" << endmsg;
    exit(1);
  }
  return 0;
}


float* Event::pid(unsigned int i)
{
  if (particles.size() > i)
    return particles[i]->pidVector;
  else {
    Alert << "accessing pid of " << i << "th of " 
	  << particles.size() << " particles" << endmsg;
    exit(1);
  }
}


void Event::addParticle(double e, double px, double py, double pz)
{
  EvtPartData* evtData = new EvtPartData(e,px,py,pz);
  particles.push_back(evtData);
  return;
}

int Event::size()
{
  return particles.size();
}
