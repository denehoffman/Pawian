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
#include "ErrLogger/ErrLogger.hh"

CBElsaReader::CBElsaReader()
{}

CBElsaReader::CBElsaReader(const std::vector<std::string>& files, int particles, int skip, bool useWeight):
  EventReader(useWeight), 
  numParticles(particles),
  linesToSkip(skip)
{
  if (0 == files.size()) {
    Alert << "empty list of event files" << endmsg;
    exit(1);
  }
  std::vector<std::string>::const_iterator iter = files.begin();
  for (; iter != files.end(); ++iter)
    fileNames.push_back(*iter);
  currentFile = fileNames.begin();
}

CBElsaReader::~CBElsaReader()
{}

bool CBElsaReader::fillAll(EventList& evtList)
{
  
  while (currentFile != fileNames.end()) {
    currentStream.open(currentFile->c_str());
    if (!currentStream) {
      Alert << "can not open " << *currentFile << endmsg;
      exit(1);
    }

    while (!currentStream.eof()) {
      double e,px,py,pz;
      Event* newEvent = new Event();
      int parts;

      if(_useWeight)
      {
         double weight;
         currentStream >> weight;
         newEvent->addWeight(weight);
      }

      for (parts = 0; parts < numParticles; parts++) {
	currentStream >> e >> px >> py >> pz;
	newEvent->addParticle(e*1.e-3,px*1.e-3,py*1.e-3,pz*1.e-3);
      }

      if (!currentStream.fail()) {
	evtList.add(newEvent);
	for (parts = 0; parts < linesToSkip; parts++)
	  currentStream >> e >> px >> py >> pz;
      }
    }
    currentStream.close();
    ++currentFile;
  }
  evtList.rewind();

  return true; // success
}



