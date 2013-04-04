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

#include "Event/CBElsaReader.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/PdtParser.hh"
#include "ErrLogger/ErrLogger.hh"

#include "qft++/topincludes/tensor.hh"
#include <vector>
#include <string>

int main()
{
  ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
  std::vector<std::string> fileNames;

  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR"); 

  ParticleTable pTable;
  PdtParser parser;
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");
  if (!parser.parse(pdtFile, pTable)) {
    Alert << "Error: could not parse " << pdtFile << endmsg;
    exit(1);
  }

  std::string pipiomegaFile(theSourcePath+"/Event/710_1350.dat");
  fileNames.push_back(pipiomegaFile);

  CBElsaReader eventReader(fileNames, 3, 1);
  EventList pipiomegaEvents;
  eventReader.fillAll(pipiomegaEvents);

  if (!pipiomegaEvents.findParticleTypes(pTable))
    Warning << "could not find all particles" << endmsg;

  Info << "\nFile has " << pipiomegaEvents.size() << " events. Each event has "
       <<  pipiomegaEvents.nextEvent()->size() << " final state particles.\n" << endmsg;
  pipiomegaEvents.rewind();

  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = pipiomegaEvents.nextEvent()) != 0 && evtCount < 20) {
    Info << "\n" 
	 << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
	 << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
	 << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
	 << endmsg;
    ++evtCount;
  }

  exit(0);
}
