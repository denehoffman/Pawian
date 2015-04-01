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

// AppBase class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <memory>

#include <cassert>
#include <memory>

#include "Event/EventList.hh"
#include "MinFunctions/AbsFcn.hh"
#include "PwaUtils/AbsChannelEnv.hh"

class fitParCol;
class NetworkClient;
class EvtDataBaseList;
class AbsPawianParameters;

class AppBase{

public:

  AppBase();
  /** Destructor */
  virtual ~AppBase();

  virtual void dumpDefaultParams();
  virtual void generate(fitParCol& theParams);
  virtual void readEvents(EventList& theEventList, std::vector<std::string>& fileNames, ChannelID channelID, bool withEvtWeight=false, int evtStart=0, int evtStop=1000000);
  virtual void createLhObjects();
  virtual void qaMode(fitParCol& theStartParams, double evtWeightSumData, int noOfFreeFitParams);
  virtual void qaModeSimple(EventList& dataEventList, EventList& mcEventList, fitParCol& theStartParams, std::shared_ptr<EvtDataBaseList> evtDataBaseList, int noOfFreeFitParams);
  virtual void plotMode(EventList& dataEventList, EventList& mcEventList, std::shared_ptr<EvtDataBaseList> evtDataBaseList);
  virtual void fixParams(std::shared_ptr<AbsPawianParameters> upar, std::vector<std::string> fixedParams);
  virtual void fixAllReleaseScaleParams(std::shared_ptr<AbsPawianParameters> upar);
  virtual bool calcAndSendClientLh(NetworkClient& theClient, fitParCol& theStartparams, ChannelID channelID);
  virtual void fitServerMode(std::shared_ptr<AbsPawianParameters> upar);
  virtual void fitNonServerMode(std::shared_ptr<AbsPawianParameters> upar, double evtWeightSumData);
  virtual void fitClientMode(fitParCol& theStartparams);
};



