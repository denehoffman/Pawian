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
// Copyright 2013 Bertram Kopf

#pragma once

#include "Event/EventReader.hh"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

class EventList;

class EventReaderDefault : public EventReader {
public:
  EventReaderDefault();
  EventReaderDefault(const std::vector<std::string> &files, int particles,
                     int skip, bool useWeight = false);

  virtual ~EventReaderDefault();

  virtual bool fill(EventList &evtList, int evtStart = 0,
                    int evtStop = 1000000);
  virtual void setUnit(const std::string &theUnit);
  virtual void setOrder(const std::string &theOrder);

protected:
  std::vector<std::string> fileNames;
  std::vector<std::string>::const_iterator currentFile;
  std::ifstream currentStream;
  int numParticles;
  int linesToSkip;
  std::string _unit;
  std::string _order;
  double _unitScaleFactor;
  bool _energyFirst;
};
