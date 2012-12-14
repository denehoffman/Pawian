#pragma once

#include "PwaUtils/EvtDataBaseListNew.hh"


class EventList;

class pbarpDataBaseList : public EvtDataBaseListNew {

public:
  pbarpDataBaseList();
  virtual ~pbarpDataBaseList();

protected:
  virtual void read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList, double& evtWeightSum, int maxEvts );

};
