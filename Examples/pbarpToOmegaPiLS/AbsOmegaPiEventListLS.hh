#ifndef _AbsOmegaPiEventListLS_H
#define _AbsOmegaPiEventListLS_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/pbarpToOmegaPiLS/OmegaPiDataLS.hh"

using OmegaPiDataLS::OmPiEvtDataLS;

class EventList;

class AbsOmegaPiEventListLS {

public:

  // create/copy/destroy:

  ///Constructor 
  AbsOmegaPiEventListLS(EventList& evtListData, EventList& evtListMc, unsigned jmax, unsigned pbarmom);


  /** Destructor */
  virtual ~AbsOmegaPiEventListLS();

  // Getters:
  const std::vector<OmPiEvtDataLS*> getDataVecs() const {return _dataList;}
  const std::vector<OmPiEvtDataLS*> getMcVecs() const {return _mcList;}
  const unsigned lMax() const {return _jmax;}
  const unsigned pbarMom() const {return _pbarmom;}

protected:
  virtual void read4Vecs(EventList& evtList, std::vector<OmPiEvtDataLS*>& omPiEvtList)=0;
  unsigned _jmax;
  unsigned _pbarmom;
  std::vector<OmPiEvtDataLS*> _dataList;
  std::vector<OmPiEvtDataLS*> _mcList;

private:


};

#endif
