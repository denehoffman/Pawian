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

#include "ErrLogger/ErrLogger.hh"
#include "Event/HepMCEventList.hh"
#include "HepMC/GenEvent.h"
#include <vector>
#include <string>

int main()
{
  ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
  std::vector<std::string> fileNames;

  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR"); 

  std::string hepMCFile(theSourcePath+"/Event/HepMCEvt.out");

  HepMCEventList hepMCEvtList(hepMCFile);
  hepMCEvtList.rewind();

  HepMC::GenEvent* genEvent;
  int evtCount = 0;
  while ( (genEvent = hepMCEvtList.nextEvent()) !=0 ){
    genEvent->print(std::cout);
    evtCount++;
  }

  Info << "HepMCEventList contains " << evtCount << " events" << endmsg;


  exit(0);
}
