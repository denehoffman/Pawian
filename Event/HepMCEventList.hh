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
#include <string>
#include "HepMC/GenEvent.h"

class HepMCEventList
{
public:
  HepMCEventList();  
  HepMCEventList(const std::vector<std::string>& files);
  HepMCEventList(const std::string& file);
  ~HepMCEventList();

  void add(HepMC::GenEvent*);
  void removeEvents(unsigned int nBegin, unsigned int nEnd);

  HepMC::GenEvent* nextEvent();
  void rewind();
  int size();

private:
  bool fillFromFiles(const std::vector<std::string>& files);
  bool fillFromFile(const std::string& file);

  std::vector<HepMC::GenEvent*> eventList;
  std::vector<HepMC::GenEvent*>::const_iterator currentEvent;
  
};
