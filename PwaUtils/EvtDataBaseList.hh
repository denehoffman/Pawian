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
  void ratioMcToData (double mcToDataRatio) {_mcToDataRatio=mcToDataRatio;}

protected:
  std::vector<EvtData*> _evtDataList;
  std::vector<EvtData*> _mcDataList;
  virtual void read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList, double& evtWeightSum, int maxEvts );
  virtual std::string getName(std::vector<Particle*>& theVec);

  double _noOfWeightedDataEvts;
  double _noOfWeightedMcEvts;
  int _mcToDataRatio;
  bool _alreadyRead;
  int _evtNoAll;
  AbsEnv* _absEnv;

};
