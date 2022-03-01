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
#include <memory>
#include <boost/unordered_map.hpp>
#include <thread>

#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsDynamics.hh"
#include "FitParams/AbsParamHandler.hh"
#include "Utils/FunctionUtils.hh"

class AbsDecay;
class FsParticleProjections;
//class AbsDynamics;
//typedef CacheVector<std::map<unsigned short, std::map<Id2StringType, complex<float> > > >  intUShortId2StringComplFloatMap;

typedef CacheVector<std::map<unsigned short, std::map<short, complex<float> > > >  intUShortShortComplFloatMap;
class AbsXdecAmp : public AbsParamHandler{

public:
  AbsXdecAmp(std::shared_ptr<AbsDecay> theDec, ChannelID channelID);
  virtual ~AbsXdecAmp();

  virtual complex<double> XdecAmp(const Spin& lamX, EvtData* theData, AbsXdecAmp* grandmaAmp=0)=0;
  virtual complex<double> XdecPartAmp(const Spin& lamX, Spin& lamDec, short fixDaughterNr,
                                      EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp=0)=0;

  virtual complex<double> daughterAmp(const Spin& lam1, const Spin& lam2, EvtData* theData);
  virtual void print(std::ostream& os) const=0;
  const std::string name() const {return _name;}
  virtual std::shared_ptr<const jpcRes>& jpcPtr() {return _JPCPtr;}
  std::shared_ptr<AbsDecay> absDec() {return _decay;}
  const std::string jpcDecsName() const {return _jpcDecsName;}
  const std::string isoKey() const {return _isoKey;}
  const std::string wignerDKey() const {return _wignerDKey;}
  //  const std::string refKey() const {return _refKey;}
  virtual void cacheAmplitudes();
  virtual bool checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld);
  virtual void calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp=0);
  virtual void setSpinProjections(const std::vector<Spin>& currentProjection, int projId);
  std::shared_ptr<AbsDynamics> getDyn() {return _absDyn;}
  std::shared_ptr<AbsXdecAmp> absXDecAmpDaughter1() {return _decAmpDaughter1;}
  std::shared_ptr<AbsXdecAmp> absXDecAmpDaughter2() {return _decAmpDaughter2;}
  //  virtual void retrieveWignerDs(EvtData* theData) { return;}

protected:

  ChannelID _channelID;
  std::shared_ptr<AbsDecay> _decay;
  const std::string _name;
  std::shared_ptr<const jpcRes> _JPCPtr;
  std::string _jpcDecsName;
  const std::string _isoKey;
  std::shared_ptr<AbsDynamics> _absDyn;
  const std::vector<std::string> _hypVec;
  std::shared_ptr<AbsXdecAmp> _decAmpDaughter1;
  std::shared_ptr<AbsXdecAmp> _decAmpDaughter2;
  Spin _Jdaughter1;
  Spin _Jdaughter2;
  std::shared_ptr<FsParticleProjections> _fsParticleProjections; 
  std::string _daughter1Name;
  std::string _daughter2Name;
  unsigned int _daughter1ProjId;
  unsigned int _daughter2ProjId;
  //  std::map<std::thread::id, unsigned short> _projIdThreadMap;
  //  std::map<std::thread::id, Spin> _lam1MinThreadMap;
  //  std::map<std::thread::id, Spin> _lam1MaxThreadMap;
  //  std::map<std::thread::id, Spin> _lam2MinThreadMap;
  //  std::map<std::thread::id, Spin> _lam2MaxThreadMap;

  unsigned short _projId;
  Spin _lam1MinProj;
  Spin _lam1MaxProj;
  Spin _lam2MinProj;
  Spin _lam2MaxProj;
  
  bool _isWeakDecay;
  double _isospinCG;
  const double _preFactor;
  std::string _key;
  const std::string _wignerDKey;
  //  std::string _refKey;

  bool _daughter1IsStable;
  bool _daughter2IsStable;
  bool _enabledlamFsDaughter1;
  bool _enabledlamFsDaughter2;

  Spin _J;
  Spin _lam1Min;
  Spin _lam1Max;
  Spin _lam2Min;
  Spin _lam2Max;

  bool _isFormationAmp;

  //  intUShortId2StringComplFloatMap _cachedAmpIdMap;
  intUShortShortComplFloatMap _cachedAmpIdMap;
  //  std::map<std::thread::id, std::map<unsigned short, complex<float> > > _cachedDynIdMap;
  std::map<unsigned short, complex<float> > _cachedDynIdMap;
  virtual void initialize();
};
