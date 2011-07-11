#ifndef _OmegaPiEventListLS_H
#define _OmegaPiEventListLS_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/pbarpToOmegaPiLS/OmegaPiDataLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"

using OmegaPiDataLS::OmPiEvtDataLS;

class EventList;

class OmegaPiEventListLS : public AbsOmegaPiEventListLS{

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaPiEventListLS(EventList& evtListData, EventList& evtListMc, unsigned jmax, unsigned pbarmom);


  /** Destructor */
  virtual ~OmegaPiEventListLS();

  // Getters:

protected:

 virtual void read4Vecs(EventList& evtList, std::vector<OmPiEvtDataLS*>& omPiEvtList);

private:

  
};

#endif
