#ifndef _JpsiGamKsKlKKEventList_H
#define _JpsiGamKsKlKKEventList_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/EvtDataBaseList.hh"



class EventList;

class JpsiGamKsKlKKEventList : public EvtDataBaseList {

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamKsKlKKEventList(EventList& evtListData, EventList& evtListMc);



  /** Destructor */
  virtual ~JpsiGamKsKlKKEventList();

  // Getters:
  

protected:


  void read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList);
  
private:


};

#endif 
