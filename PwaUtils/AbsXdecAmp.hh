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

// AbsXdecAmp class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <mutex>
#include <memory>
#include <boost/unordered_map.hpp>

#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/AbsParamHandler.hh"
#include "PwaUtils/AbsDynamics.hh"

class AbsDecay;
//class AbsDynamics;

typedef boost::unordered_map<int, std::map< unsigned short, complex<float> > > intShortComplFloatMap;
//typedef boost::unordered_map<int, std::map<Spin, std::map<Spin, complex<float> > > > intSpinSpinFloatUsMap;


class AbsXdecAmp : public AbsParamHandler{

public:
  AbsXdecAmp(std::shared_ptr<AbsDecay> theDec, ChannelID channelID);
  virtual ~AbsXdecAmp();

  virtual complex<double> XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp=0)=0;
  virtual complex<double> XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr,
                                      EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp=0)=0;

  virtual complex<double> daughterAmp(Spin& lam1, Spin& lam2, EvtData* theData, Spin& lamFs);
  virtual void print(std::ostream& os) const=0;
  const std::string name() const {return _name;}
  virtual std::shared_ptr<const jpcRes>& jpcPtr() {return _JPCPtr;}
  std::shared_ptr<AbsDecay> absDec() {return _decay;}
  const std::string jpcDecsName() const {return _jpcDecsName;}
  const std::string isoKey() const {return _isoKey;}

  virtual void cacheAmplitudes();

protected:

  ChannelID _channelID;
  std::shared_ptr<AbsDecay> _decay;
  const std::string _name;
  std::shared_ptr<const jpcRes> _JPCPtr;
  const std::string _jpcDecsName;
  const std::string _isoKey;
  std::shared_ptr<AbsDynamics> _absDyn;
  const std::vector<std::string> _hypVec;
  std::shared_ptr<AbsXdecAmp> _decAmpDaughter1;
  std::shared_ptr<AbsXdecAmp> _decAmpDaughter2;
  Spin _Jdaughter1;
  Spin _Jdaughter2;
  double _isospinCG;
  const double _preFactor;
  std::string _key;
  const std::string _wignerDKey;
  std::mutex theMutex;

  bool _daughter1IsStable;
  bool _daughter2IsStable;
  bool _enabledlamFsDaughter1;
  bool _enabledlamFsDaughter2;

  intShortComplFloatMap _cachedAmpShortMap;
  //  intSpinSpinFloatUsMap _cachedAmpMap;

  virtual void initialize();
};
