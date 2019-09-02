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

// EvtDataScatteringList class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#pragma once

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <memory>
#include <boost/unordered_map.hpp>
#include <map>

#include "PwaUtils/EvtDataBaseList.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "Particle/Particle.hh"

class EventList;
class Event;
class Particle;

class EvtDataScatteringList : public EvtDataBaseList{

public:
  EvtDataScatteringList(ChannelID channelID);
  virtual ~EvtDataScatteringList();

  virtual void read(EventList& evtListData, EventList& evtListMc);

  virtual void readScatteringData(EventList& evtList, std::vector<EvtData*>& theEvtList, int startNo);
  virtual void readScatteringDataInterpol(EventList& evtList, std::vector<EvtData*>& theEvtList, int startNo);
  virtual void readScatteringDefaultFit(std::vector<EvtData*>& dataEvtList, std::vector<EvtData*>& fitEvtList, int startNo);
  virtual EvtData* convertEvent(Event* theEvent, int evtNo=1);
  virtual EvtData* convertEventInterpol(Event* oldEvent, Event* newEvent, int evtNo);
  virtual std::vector<EvtData*> convertEventInterpol(Event* oldEvent, Event* newEvent, int evtNo, unsigned int noOfPoints);
 
  static std::string M_PIPISCAT_NAME;
  static std::string DATA_PIPISCAT_NAME;
  static std::string DATAERR_PIPISCAT_NAME;
  static std::string FIT_PIPISCAT_NAME;
  static std::string FITERR_PIPISCAT_NAME;

protected:
  int _noOfInterPolData;
};
