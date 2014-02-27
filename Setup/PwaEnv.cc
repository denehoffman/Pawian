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
#include "Particle/PdtParser.hh"
#include "Setup/SetupParser.hh"
#include "Event/EventList.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Event/CBElsaReader.hh"
#include "DecayTree/DecayTree.hh"

PwaEnv& PwaEnv::instance()
{
  static PwaEnv theEnv;
  return theEnv;
}

bool PwaEnv::setup(std::string& setupFileName)
{
  PdtParser pdtParser;
  std::string theSourcePath=getenv("TOP_DIR"); 
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");
  theParticleTable = new ParticleTable;

  if (!pdtParser.parse(pdtFile, *theParticleTable)) {
    Alert << "can not parse particle table " << pdtFile << endmsg;
    exit(1);
  }

  SetupParser setupParser;
  if (!setupParser.parse(setupFileName, theParticleTable)) {
    Alert << "can not parse setup " << setupFileName << endmsg;
    exit(1);
  }

  theDecayTree = new DecayTree(setupParser.edgeList());
  if (0 == theDecayTree) {
    Alert << "can not build decay tree" << endmsg;
    exit(1);
  }

  if (setupParser.setup()->beamInput.size() == 0) {
    Alert << "no beam data defined in setup" << endmsg;
    exit(1);
  }

  CBElsaReader eventReader(setupParser.setup()->beamInput, 3, 1);
  theBeamEventList = new EventList();
  eventReader.fillAll(*theBeamEventList);

  if (setupParser.setup()->mcInput.size() == 0) {
    Alert << "no MC data defined in setup" << endmsg;
    exit(1);
  }

  CBElsaReader mcEventReader(setupParser.setup()->mcInput, 3, 1);
  theMcEventList = new EventList();
  mcEventReader.fillAll(*theMcEventList);

  return 1;
}

ParticleTable* PwaEnv::particleTable()
{
  return theParticleTable;
}

EventList* PwaEnv::beamEventList()
{
  return theBeamEventList;
}

EventList* PwaEnv::mcEventList()
{
  return theMcEventList;
}

PwaEnv::PwaEnv()
{
}

PwaEnv::~PwaEnv()
{
  delete theParticleTable;
  delete theBeamEventList;
  delete theMcEventList;
}
