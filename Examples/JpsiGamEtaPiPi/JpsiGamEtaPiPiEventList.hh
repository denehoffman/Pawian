#ifndef _JpsiGamEtaPiPiEventList_H
#define _JpsiGamEtaPiPiEventList_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/EvtDataBaseList.hh"



class EventList;

class JpsiGamEtaPiPiEventList : public EvtDataBaseList {

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamEtaPiPiEventList(EventList& evtListData, EventList& evtListMc);



  /** Destructor */
  virtual ~JpsiGamEtaPiPiEventList();

  // Getters:
  

protected:


  void read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList);
  
private:


};

#endif 
