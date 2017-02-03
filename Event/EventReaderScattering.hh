//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#pragma once

#include "Event/EventReaderDefault.hh"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <utility>


class EventList;

class EventReaderScattering : public EventReaderDefault
{
public:
  EventReaderScattering();
  EventReaderScattering(const std::vector<std::string>& files, int particles, int skip, bool useWeight=false);

  virtual ~EventReaderScattering();

  virtual bool fill(EventList& evtList, int evtStart=0, int evtStop=1000000);
  virtual void setOrder(const std::string& theOrder);
  
protected:
};


