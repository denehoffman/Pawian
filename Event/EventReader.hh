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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Event/MassRangeCut.hh"

class EventList;

class EventReader
{
public:
  EventReader();
  EventReader(bool useWeight);
  virtual ~EventReader();

  virtual bool fill(EventList& evtList, int evtStart=0, int evtStop=1000000) = 0;

  virtual void setMassRange(std::vector< std::shared_ptr<MassRangeCut> >& massRangeCuts){
    _useMassRange=true;
    _massRangeCuts=massRangeCuts;
  };


protected:
  bool _useWeight;
  bool _useMassRange;
  std::vector< std::shared_ptr<MassRangeCut> > _massRangeCuts;
};
