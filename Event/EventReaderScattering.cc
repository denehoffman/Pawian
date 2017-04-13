//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// EventReaderScattering class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include "Event/EventReaderScattering.hh"

#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"

EventReaderScattering::EventReaderScattering() :
  EventReaderDefault()
{
  _order=std::string("Mass Phi PhiErr Eta EtaErr");
  numParticles=1;
}

EventReaderScattering::EventReaderScattering(const std::vector<std::string>& files, int particles, int skip, bool useWeight):
  EventReaderDefault(files, particles, skip, useWeight)
{
  _order=std::string("Mass Phi PhiErr Eta EtaErr");
  numParticles=1;
}

EventReaderScattering::~EventReaderScattering()
{}

bool EventReaderScattering::fill(EventList& evtList, int evtStart, int evtStop)
{
  int currentEvtNo=-1;

  while (currentFile != fileNames.end()) {
    currentStream.open(currentFile->c_str());
    if (!currentStream) {
      Alert << "can not open " << *currentFile ;  // << endmsg;
      exit(1);
    }

    while (!currentStream.eof()) {
      currentEvtNo++;
      // double m,phi,phiErr,mu,muErr;
      double m, dataPoint, dataPointErr;      
      Event* newEvent = new Event(currentEvtNo);
      
      currentStream >> m >> dataPoint >> dataPointErr;
      double scaledMass=m/_unitScaleFactor;
      newEvent->addParticle(scaledMass, 0., 0., 0.);
      newEvent->addScatterInfo(dataPoint, dataPointErr);
      Vector4<double> tmp = newEvent->p4(0);

      bool acceptEvt=true;
      if(_useMassRange){
	std::vector< std::shared_ptr<MassRangeCut> >::iterator itMassRangeCut;
	for (itMassRangeCut=_massRangeCuts.begin(); itMassRangeCut!=_massRangeCuts.end(); ++itMassRangeCut){
	  if(scaledMass < (*itMassRangeCut)->massMin() || scaledMass > (*itMassRangeCut)->massMax()){
	    acceptEvt=false;
	    break;
	  } 
	}
      }

      if( !acceptEvt || currentEvtNo<evtStart || currentEvtNo>evtStop){
	if (!acceptEvt) currentEvtNo--;
	delete newEvent;
	continue;
      }
 
      if (!currentStream.fail()) {
	InfoMsg << "added data point (event) for pi pi scattering no: " << newEvent->eventNo() 
		<< "\nm: " << newEvent->p4(0)->Mass() 
		<< "\ndataPoint: " << newEvent->DataPoint() <<  "\t error: " << newEvent->DataPointErr()
		<< "\n" << endmsg; 
      	evtList.add(newEvent);
      	// for (parts = 0; parts < linesToSkip; parts++)
      	//   currentStream >> m >> phi >> phiErr >> mu >> muErr;
	
      }
    }
    if (currentEvtNo == 0){
      Alert << "No events have been read in from File: " << *currentFile << " ! \nNo events in massRange? Wrong order of particles? currentEvtNo = " << currentEvtNo;
      exit(1);
    }
    currentStream.close();
    ++currentFile;
  }
  evtList.rewind();
 
  return true; // success
}

void EventReaderScattering::setOrder(const std::string& theOrder){
}

