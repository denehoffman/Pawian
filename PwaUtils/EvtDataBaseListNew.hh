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

struct EvtDataNew {
  mapInt4Vec FourVecsProd;
  mapInt4Vec FourVecsDec;
  mapIntSpinComplex WignerDsProd;
  mapIntSpinComplex WignerDsDec;
  double evtWeight;
};

class EventList;

class EvtDataBaseListNew {

public:
  EvtDataBaseListNew();
  virtual ~EvtDataBaseListNew();

  void read(EventList& evtListData, EventList& evtListMc);

  const std::vector<EvtDataNew*> getDataVecs() const { return _evtDataList; }
  const std::vector<EvtDataNew*> getMcVecs() const { return _mcDataList; }

  double NoOfWeightedDataEvts() const {return _noOfWeightedDataEvts;}
  double NoOfWeightedMcEvts() const {return _noOfWeightedMcEvts;}
  void ratioMcToData (double mcToDataRatio) {_mcToDataRatio=mcToDataRatio;}

protected:
  std::vector<EvtDataNew*> _evtDataList;
  std::vector<EvtDataNew*> _mcDataList;
  virtual void read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList, double& evtWeightSum, int maxEvts ) = 0;

  double _noOfWeightedDataEvts;
  double _noOfWeightedMcEvts;
  int _mcToDataRatio;
  bool _alreadyRead;

};
