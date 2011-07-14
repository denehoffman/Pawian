#ifndef _OmegaTo3PiEventListLS_H
#define _OmegaTo3PiEventListLS_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/pbarpToOmegaPiLS/OmegaPiDataLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"

using OmegaPiDataLS::OmPiEvtDataLS;

class EventList;

class OmegaTo3PiEventListLS : public AbsOmegaPiEventListLS{

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaTo3PiEventListLS(EventList& evtListData, EventList& evtListMc, unsigned jmax, unsigned pbarmom);


  /** Destructor */
  virtual ~OmegaTo3PiEventListLS();

  // Getters:

protected:

 virtual void read4Vecs(EventList& evtList, std::vector<OmPiEvtDataLS*>& omPiEvtList);

private:

  
};

#endif
