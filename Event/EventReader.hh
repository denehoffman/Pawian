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

class EventList;

class EventReader
{
public:
  EventReader();
  EventReader(bool useWeight);
  virtual ~EventReader();

  virtual bool fill(EventList& evtList, int evtStart=0, int evtStop=1000000) = 0;

  virtual void setMassRange(EventList& evtList, double massMin, double massMax, std::vector<unsigned int>& particleIndices){
    _useMassRange=true;
    _massMin = massMin;
    _massMax = massMax;
    _particleIdxesMassRange=particleIndices;
    //   return fillAll(evtList);
  };

protected:
  bool _useWeight;
  bool _useMassRange;
  double _massMin;
  double _massMax;
  std::vector<unsigned int> _particleIdxesMassRange;

  bool isMassrangeParticle(int index); 
};
