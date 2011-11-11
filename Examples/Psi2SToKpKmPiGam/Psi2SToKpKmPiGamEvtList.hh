#ifndef _Psi2SToKpKmPiGamEvtList_H
#define _Psi2SToKpKmPiGamEvtList_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "PwaUtils/EvtDataBaseList.hh"

class EventList;

class Psi2SToKpKmPiGamEvtList : public EvtDataBaseList {

public:

  // create/copy/destroy:

  ///Constructor 
  Psi2SToKpKmPiGamEvtList(EventList& evtListData, EventList& evtListMc);



  /** Destructor */
  virtual ~Psi2SToKpKmPiGamEvtList();

  // Getters:


protected:
  virtual void read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList);

private:

};

#endif 
