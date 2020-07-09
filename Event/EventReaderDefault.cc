//************************************************************************//
//                                                                        //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)                      //
//               Julian Pychy (julian@ep1.rub.de)                         //
//               - Ruhr-Universität Bochum                                //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

// EventReaderDefault class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include "Event/EventReaderDefault.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Event/Event.hh"
#include "Event/EventList.hh"

EventReaderDefault::EventReaderDefault()
    : EventReader(false), numParticles(0), linesToSkip(0), _unit("GEV"),
      _order("Px Py Pz E"), _unitScaleFactor(1.), _energyFirst(false) {}

EventReaderDefault::EventReaderDefault(const std::vector<std::string> &files,
                                       int particles, int skip, bool useWeight)
    : EventReader(useWeight), numParticles(particles), linesToSkip(skip),
      _unit("GEV"), _order("Px Py Pz E"), _unitScaleFactor(1.),
      _energyFirst(false) {
  if (0 == files.size()) {
    Alert << "empty list of event files"; // << endmsg;
    exit(1);
  }
  std::vector<std::string>::const_iterator iter = files.begin();
  for (; iter != files.end(); ++iter)
    fileNames.push_back(*iter);
  currentFile = fileNames.begin();
}

EventReaderDefault::~EventReaderDefault() {}

bool EventReaderDefault::fill(EventList &evtList, int evtStart, int evtStop) {
  int currentEvtNo = -1;

  while (currentFile != fileNames.end()) {
    currentStream.open(currentFile->c_str());
    if (!currentStream.is_open()) {
      Alert << "can not open " << *currentFile << endmsg;
      exit(1);
    }

    while (!currentStream.eof()) {
      currentEvtNo++;
      double e, px, py, pz;
      Event *newEvent = new Event();
      int parts;

      if (_useWeight) {
        double weight;
        currentStream >> weight;
        newEvent->addWeight(weight);
      }
      Vector4<double> fvX(0, 0, 0, 0); // X four-vector
      for (parts = 0; parts < numParticles; parts++) {
        if (_energyFirst)
          currentStream >> e >> px >> py >> pz;
        else
          currentStream >> px >> py >> pz >> e;
        newEvent->addParticle(e / _unitScaleFactor, px / _unitScaleFactor,
                              py / _unitScaleFactor, pz / _unitScaleFactor);
        Vector4<double> tmp = newEvent->p4(parts);
      }

      bool acceptEvt = true;
      if (_useMassRange) {
        std::vector<std::shared_ptr<MassRangeCut>>::iterator itMassRangeCut;
        for (itMassRangeCut = _massRangeCuts.begin();
             itMassRangeCut != _massRangeCuts.end(); ++itMassRangeCut) {
          Vector4<float> particleSystem4Vec;
          std::vector<unsigned int> particleIndices =
              (*itMassRangeCut)->particleIds();
          for (auto it = particleIndices.begin(); it != particleIndices.end();
               ++it) {
            particleSystem4Vec += *(newEvent->p4(*it));
          }
          double invMass = particleSystem4Vec.Mass();
          if (!((*itMassRangeCut)->isAntiCut())) {
            if (invMass < (*itMassRangeCut)->massMin() ||
                invMass > (*itMassRangeCut)->massMax()) {
              acceptEvt = false;
              break;
            }
          } else { // anti-cut
            if (invMass > (*itMassRangeCut)->massMin() &&
                invMass < (*itMassRangeCut)->massMax()) {
              acceptEvt = false;
              break;
            }
          }
        }
      }

      if (!acceptEvt || currentEvtNo < evtStart || currentEvtNo > evtStop) {
        if (!acceptEvt)
          currentEvtNo--;
        delete newEvent;
        continue;
      }

      if (!currentStream.fail()) {
        evtList.add(newEvent);
        for (parts = 0; parts < linesToSkip; parts++)
          if (_energyFirst)
            currentStream >> e >> px >> py >> pz;
          else
            currentStream >> px >> py >> pz >> e;
      }
    }
    if (currentEvtNo < 0) {
      Alert << "No events have been read in from File: " << *currentFile
            << " ! \nNo events in massRange? Wrong order of particles? "
               "currentEvtNo = "
            << currentEvtNo;
      exit(1);
    }
    currentStream.close();
    ++currentFile;
  }
  evtList.rewind();

  return true; // success
}

void EventReaderDefault::setUnit(const std::string &theUnit) {
  _unit = theUnit;
  if (_unit == "GEV") {
    _unitScaleFactor = 1.;
  } else if (_unit == "MEV") {
    _unitScaleFactor = 1000.;
  } else {
    Alert << "unit " << _unit << " does not exist!!!" << endmsg;
    exit(0);
  }
}

void EventReaderDefault::setOrder(const std::string &theOrder) {
  _order = theOrder;
  if (_order == "E Px Py Pz") {
    _energyFirst = true;
  } else if (_order == "Px Py Pz E") {
    _energyFirst = false;
  } else {
    Alert << "order " << _order << " does not exist!!!" << endmsg;
    exit(0);
  }
}
