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

#include "Event/EventList.hh"
#include "PwaUtils/AbsFcn.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"

// #include <memory>

class AbsLh;
class fitParams;
class AbsHist;
class NetworkClient;
class EvtDataBaseList;
class AbsHist;

using namespace ROOT::Minuit2;

class AppBase{

public:

  AppBase();
  /** Destructor */
  virtual ~AppBase();

  virtual void dumpDefaultParams();
  virtual void generate(fitParams& theParams);
  virtual void readEvents(EventList& theEventList, std::vector<std::string>& fileNames, ChannelID channelID, bool withEvtWeight=false, int evtStart=0, int evtStop=1000000);
  virtual void qaMode(fitParams& theStartParams, double evtWeightSumData, int noOfFreeFitParams);
  virtual void qaModeSimple(EventList& dataEventList, EventList& mcEventList, fitParams& theStartParams, std::shared_ptr<EvtDataBaseList> evtDataBaseList, int noOfFreeFitParams);
  virtual void plotMode(EventList& dataEventList, EventList& mcEventList, std::shared_ptr<EvtDataBaseList> evtDataBaseList);
  virtual void fixParams(MnUserParameters& upar, std::vector<std::string> fixedParams);
  virtual void fixAllReleaseScaleParams(MnUserParameters& upar);
  virtual FunctionMinimum migradDefault(AbsFcn& theFcn, MnUserParameters& upar);
  virtual void printFitResult(FunctionMinimum& min, fitParams& theStartparams, std::ostream& os, double evtWeightSumData=0, int noOfFreeFitParams=0);
  virtual bool calcAndSendClientLh(NetworkClient& theClient, fitParams& theStartparams, ChannelID channelID);

};



