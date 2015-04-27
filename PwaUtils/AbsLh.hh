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

// AbsLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <memory>

#include "PwaUtils/AbsParamHandler.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "FitParams/FitParColBase.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/DataUtils.hh"


class AbsLh : public AbsParamHandler{

public:
  AbsLh(std::shared_ptr<AbsLh>);
  AbsLh(ChannelID channelID);
  void initialize();

  virtual ~AbsLh();
  virtual AbsLh* clone_() const = 0;

  virtual double calcLogLh(fitParCol& theParamVal);
  static  double mergeLogLhData(LHData& theLHData);
  virtual double addDataToLogLh(EvtData* dataEvt, fitParCol& theParamVal, LHData& theLHData);
  virtual double addMcToLogLh(EvtData* mcEvt, fitParCol& theParamVal, LHData& theLHData);
  virtual void calcLogLhDataClient(fitParCol& theParamVal, LHData& theLHData);
  virtual double calcEvtIntensity(EvtData* theData, fitParCol& theParamVal)=0;

  virtual void setDataVec(std::vector<EvtData*> theVec);
  virtual void setMcVec(std::vector<EvtData*> theVec);

  virtual std::vector<EvtData*> getDataVec() {return _evtDataVec;}
  virtual std::vector<EvtData*> getMcVec() {return _evtMCVec;}

  virtual void getDefaultParams(fitParCol& fitVal, fitParCol& fitErr);
  
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  
  virtual bool checkRecalculation(fitParCol& theParamVal);
  virtual void cacheAmplitudes();
  virtual void updateFitParams(fitParCol& theParamVal);

  virtual void print(std::ostream& os) const=0;

  ChannelID getChannelID(){return _channelID;}
  std::string getChannelScaleParam(){return _channelScaleParam;}

protected:

  ChannelID _channelID;
  std::string _channelScaleParam;
  std::vector<EvtData*> _evtDataVec;
  std::vector<EvtData*> _evtMCVec;
  std::vector< std::shared_ptr<AbsXdecAmp> > _decAmps;

  bool _usePhasespace;
  const std::string _phasespaceKey;
  std::map<const std::string, bool> _hypMap;
  unsigned int _calcCounter;
  unsigned short _noOfThreads;

  virtual void setHyps( const std::map<const std::string, bool>& theMap,
			bool& theHyp, std::string& theKey);

  virtual void ThreadfuncData(unsigned int minEvent, unsigned int maxEvent,
			      LHData& theLHData, fitParCol& theParamVal);
  virtual void ThreadfuncMc(unsigned int minEvent, unsigned int maxEvent,
			     LHData& theLHData, fitParCol& theParamVal);
};
