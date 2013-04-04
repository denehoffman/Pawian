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

#include "Setup/PwaEnv.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "qft++/topincludes/tensor.hh"

#include "ErrLogger/ErrLogger.hh"

#include <iostream>
#include <string>
#include <cstdlib>

int main()
{
  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR"); 
  std::string setupFile(theSourcePath + "/Setup/test.setup");
  PwaEnv::instance().setup(setupFile);

  ParticleTable* pTable = PwaEnv::instance().particleTable();
  if (0 == pTable) {
    Alert << "getting ParticleTable failed" << endmsg;
    exit(1);
  }
  pTable->print(std::cout);

  EventList* eventList = PwaEnv::instance().beamEventList();
  if (0 == eventList) {
    Alert << "getting beam EventList failed" << endmsg;
    exit(1);
  }

  Info << "Input file has " << eventList->size() << " events. Each event has "
       <<  eventList->nextEvent()->size() << " final state particles.\n" << endmsg;
  eventList->rewind();

  EventList* mcEventList = PwaEnv::instance().mcEventList();
  if (0 == mcEventList) {
    Alert << "getting MC EventList failed" << endmsg;
    exit(1);
  }

  Info << "MC Input file has " << mcEventList->size() << " events. Each event has "
       <<  mcEventList->nextEvent()->size() << " final state particles.\n" << endmsg;
  mcEventList->rewind();

  Event* anEvent;
  int evtCount = 0;
  Info << "======== beam events ========" << endmsg;
  while ((anEvent = eventList->nextEvent()) != 0 && evtCount < 20) {
    Info << "\n" 
	 << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
	 << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
	 << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
	 << endmsg;
    ++evtCount;
  }

  evtCount = 0;
  Info << "======== MC events ========" << endmsg;
  while ((anEvent = mcEventList->nextEvent()) != 0 && evtCount < 20) {
    Info << "\n" 
	 << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
	 << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
	 << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
	 << endmsg;
    ++evtCount;
  }


  return 0;
}
