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

// EvtDataBaseList class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <memory>
#include <boost/unordered_map.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "Particle/Particle.hh"

//typedef unsigned int Id3StringType;

typedef std::map<std::string, Vector4<double> > mapString4Vec;
//typedef std::map<std::string, map<std::string, map<Spin,map<Spin,map<Spin,complex<double> > > > > > mapStringStringSpinComplex;
//typedef std::map<std::string, std::map<std::string, std::map<Id3StringType, complex<double> > > > mapStringStringIntComplex;
typedef std::map<std::string, double> mapStringDouble;
typedef std::map<std::string, map<std::string, double> > mapStringStringDouble;
typedef std::map<std::string, map<Spin,map<Spin, complex<double> > > > mapStringComplDouble;
typedef std::map<std::string, map<Spin,map<Spin,map<Spin, map<Spin, map<Spin,complex<double> > > > > > > mapString5SpinsComplex;
typedef std::map<std::string, map<int,map<Spin, map<Spin, map<Spin,complex<double> > > > > > mapStringInt3SpinsComplex;
typedef std::map<unsigned int, std::map<Id3StringType, complex<double> > > mapIntInt3Complex;
typedef std::map<unsigned int, std::map<Id1StringType, complex<double> > > mapIntInt1Complex;
typedef std::map<unsigned int, double> mapIntDouble;
struct EvtData {
  mapString4Vec FourVecsString;
  //  mapStringStringSpinComplex WignerDsStringString;
  //  mapStringStringIntComplex WignerDStringStringId;
  mapStringDouble DoubleString;
  mapStringStringDouble DoubleStringString;
  mapStringComplDouble ComplexDoubleString;
  mapString5SpinsComplex ComplexDouble5SpinString;
  mapStringInt3SpinsComplex ComplexDoubleInt3SpinString;
  mapIntInt3Complex WignerDIdId3;
  mapIntInt1Complex WignerDIdId1;
  mapIntDouble DoubleId;
  double evtWeight;
  int evtNo;
};

class EventList;
class Event;
class Particle;

class EvtDataBaseList {

public:
  EvtDataBaseList(ChannelID channelID);
  virtual ~EvtDataBaseList();

  void read(EventList& evtListData, EventList& evtListMc);


  const std::vector<EvtData*> getDataVecs() const { return _evtDataList; }
  const std::vector<EvtData*> getMcVecs() const { return _mcDataList; }

  double NoOfWeightedDataEvts() const {return _noOfWeightedDataEvts;}
  double NoOfWeightedMcEvts() const {return _noOfWeightedMcEvts;}
  virtual void read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList, double& evtWeightSum, int maxEvts, int startNo);
  virtual EvtData* convertEvent(Event* theEvent, int evtNo=1);
  static double noOfWeightedEvts(EventList& evtList, ChannelID channelID, int maxEvts, int startNo);
 
protected:
  ChannelID _channelID;
  std::vector<EvtData*> _evtDataList;
  std::vector<EvtData*> _mcDataList;
  virtual std::string getName(std::vector<Particle*>& theVec);

  double _noOfWeightedDataEvts;
  double _noOfWeightedMcEvts;
  //  int _mcToDataRatio;
  bool _alreadyRead;
  std::vector<Particle*>  _finalStateParticles;
 
};
