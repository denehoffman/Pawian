//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

#include "Event/EventList.hh"

#include "Event/Event.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "ErrLogger/ErrLogger.hh"
#include "qft++/topincludes/tensor.hh"

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
  else {
    Alert << "can not add 0 pointer to event list" << endmsg;
    exit(1);
  }

  return;
}

void EventList::removeEvents(unsigned int nBegin, unsigned int nEnd)
{
  if ( nBegin > nEnd || eventList.size()<nEnd) {
    Alert << "can not remove event no " << nBegin << " - " << nEnd 
	  << " from list" << endmsg; 
    exit(1);  
  }

  eventList.erase(eventList.begin()+nBegin, eventList.begin()+nEnd);
}

void EventList::removeAndDeleteEvents(unsigned int nBegin, unsigned int nEnd)
{
  if ( nBegin > nEnd || eventList.size()<nEnd) {
    Alert << "can not remove event no " << nBegin << " - " << nEnd 
	  << " from list" << endmsg; 
    exit(1);  
  }

  std::vector<Event*>::iterator currentNonConstEvent;
  for (currentNonConstEvent = eventList.begin();
        currentNonConstEvent != eventList.end();
        ++currentNonConstEvent){
     delete *currentNonConstEvent;
     (*currentNonConstEvent)=0;
   }
  eventList.erase(eventList.begin()+nBegin, eventList.begin()+nEnd);
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

Particle* EventList::particle(unsigned int num)
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
  if (eventList.size() == 0) {
    Alert << "eventList is empty" << endmsg;
    exit(1);
  }

  bool result = true;
  int partsPerEvent = eventList[0]->size();
  int i;
  for (i=0; i<partsPerEvent; i++) {
    Particle* thisParticle;
    thisParticle = pdtTable.particle(eventList[0]->p4(i)->Mass());
    particleRefs.push_back(thisParticle);
    if (0 == thisParticle) {
      WarningMsg << "did not find a particle with mass " 
	      << eventList[0]->p4(i)->Mass() << endmsg;
      result = false;
    }
  }

  return result;
}
