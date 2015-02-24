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

#pragma once

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
  void removeEvents(unsigned int nBegin, unsigned int nEnd);
  void removeAndDeleteEvents(unsigned int nBegin, unsigned int nEnd);

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
