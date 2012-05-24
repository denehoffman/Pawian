#ifndef _JpsiGamEtaPiPiEventListNew_H
#define _JpsiGamEtaPiPiEventListNew_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"



class EventList;

class JpsiGamEtaPiPiEventListNew : public EvtDataBaseListNew {

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamEtaPiPiEventListNew(EventList& evtListData, EventList& evtListMc);



  /** Destructor */
  virtual ~JpsiGamEtaPiPiEventListNew();

  // Getters:
  

protected:


  void read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList);
  
private:


};

#endif 
