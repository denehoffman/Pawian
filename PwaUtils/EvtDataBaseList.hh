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
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"
#include "Particle/Particle.hh"

typedef std::map<int, Vector4<double> > mapInt4Vec;
typedef std::map<int, map<Spin,map<Spin,map<Spin,complex<double> > > > > mapIntSpinComplex;
typedef std::map<std::string, Vector4<double> > mapString4Vec;
typedef std::map<std::string, map<Spin,map<Spin,map<Spin,complex<double> > > > > mapStringSpinComplex;
typedef std::map<std::string, double> mapStringDouble;

struct EvtData {
  mapInt4Vec FourVecsProd;
  mapInt4Vec FourVecsDec;
  mapIntSpinComplex WignerDsProd;
  mapIntSpinComplex WignerDsDec;
  mapString4Vec FourVecsString;
  mapStringSpinComplex WignerDsString;
  mapStringDouble DoubleString;
  double evtWeight;
  int evtNo;
};

class EventList;
class AbsEnv;

class EvtDataBaseList {

public:
  EvtDataBaseList(AbsEnv* theEnv);
  virtual ~EvtDataBaseList();

  void read(EventList& evtListData, EventList& evtListMc);


  const std::vector<EvtData*> getDataVecs() const { return _evtDataList; }
  const std::vector<EvtData*> getMcVecs() const { return _mcDataList; }

  double NoOfWeightedDataEvts() const {return _noOfWeightedDataEvts;}
  double NoOfWeightedMcEvts() const {return _noOfWeightedMcEvts;}
  //  void ratioMcToData (double mcToDataRatio) {_mcToDataRatio=mcToDataRatio;}
  virtual void read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList, double& evtWeightSum, int maxEvts, int startNo);

protected:
  std::vector<EvtData*> _evtDataList;
  std::vector<EvtData*> _mcDataList;
  virtual std::string getName(std::vector<Particle*>& theVec);

  double _noOfWeightedDataEvts;
  double _noOfWeightedMcEvts;
  //  int _mcToDataRatio;
  bool _alreadyRead;
  AbsEnv* _absEnv;

};
