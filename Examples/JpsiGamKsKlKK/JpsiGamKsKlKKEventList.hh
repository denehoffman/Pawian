#ifndef _JpsiGamKsKlKKEventList_H
#define _JpsiGamKsKlKKEventList_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKData.hh"

using JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData;

class EventList;

class JpsiGamKsKlKKEventList {

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamKsKlKKEventList(EventList& evtListData, EventList& evtListMc);



  /** Destructor */
  virtual ~JpsiGamKsKlKKEventList();

  // Getters:
  const std::vector<JpsiGamKsKlKKEvtData*> getDataVecs() const {return _dataList;}
  const std::vector<JpsiGamKsKlKKEvtData*> getMcVecs() const {return _mcList;}


protected:


private:
  std::vector<JpsiGamKsKlKKEvtData*> _dataList;
  std::vector<JpsiGamKsKlKKEvtData*> _mcList;

  void read4Vecs(EventList& evtList, std::vector<JpsiGamKsKlKKEvtData*>& theEvtList);
  

};

#endif 
