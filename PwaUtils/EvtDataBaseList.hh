#pragma once

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <boost/shared_ptr.hpp>
#include <map>

#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"

struct EvtData {
  std::map<int, Vector4<double> > FourVecs;
  std::map<int, map<Spin,map<Spin,map<Spin,complex<double> > > > > WignerDs;
  double evtWeight;
};

class EventList;

class EvtDataBaseList {

public:
  EvtDataBaseList();
  virtual ~EvtDataBaseList();

  const std::vector<EvtData*> getDataVecs() const {return _evtDataList;}
  const std::vector<EvtData*> getMcVecs() const {return _mcDataList;}

protected:
  std::vector<EvtData*> _evtDataList;
  std::vector<EvtData*> _mcDataList;
  virtual void read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList)=0;  
};
