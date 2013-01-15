// EvtDataBaseList class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"

struct enumProd4V {
  enum { Psi=0, n4Vecs };
  
  static const std::string& name(unsigned int t) {
    static std::string fitName[enumProd4V::n4Vecs] = {"Psi"};
    if (t<0 || t>=enumProd4V::n4Vecs) assert(0);
    return fitName[t];
  }
};

struct enumProdDfunc {
  enum {Psi=0, nDfuncts};

  static const std::string& name(unsigned int t) {
    static std::string fitName[enumProdDfunc::nDfuncts] = {"Psi"};
    if (t<0 || t>=enumProdDfunc::nDfuncts) assert(0);
    return fitName[t];
  }
};

typedef std::map<int, Vector4<double> > mapInt4Vec;
typedef std::map<int, map<Spin,map<Spin,map<Spin,complex<double> > > > > mapIntSpinComplex;
typedef std::map<std::string, Vector4<double> > mapString4Vec;
typedef std::map<std::string, map<Spin,map<Spin,map<Spin,complex<double> > > > > mapStringSpinComplex;

struct EvtData {
  mapInt4Vec FourVecsProd;
  mapInt4Vec FourVecsDec;
  mapIntSpinComplex WignerDsProd;
  mapIntSpinComplex WignerDsDec;
  mapString4Vec FourVecsString;
  mapStringSpinComplex WignerDsString;
  double evtWeight;
  int evtNo;
};

class EventList;

class EvtDataBaseList {

public:
  EvtDataBaseList();
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
  virtual void read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList, double& evtWeightSum, int maxEvts ) = 0;

  double _noOfWeightedDataEvts;
  double _noOfWeightedMcEvts;
  int _mcToDataRatio;
  bool _alreadyRead;
  int _evtNoAll;

};
