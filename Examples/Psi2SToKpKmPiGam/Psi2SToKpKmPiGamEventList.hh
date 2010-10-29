#ifndef _Psi2SToKpKmPiGamEventList_H
#define _Psi2SToKpKmPiGamEventList_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"

using Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData;

class EventList;

class Psi2SToKpKmPiGamEventList {

public:

  // create/copy/destroy:

  ///Constructor 
  Psi2SToKpKmPiGamEventList(EventList& evtListData, EventList& evtListMc);



  /** Destructor */
  virtual ~Psi2SToKpKmPiGamEventList();

  // Getters:
  const std::vector<Psi2SToKpKmPiGamEvtData*> getDataVecs() const {return _dataList;}
  const std::vector<Psi2SToKpKmPiGamEvtData*> getMcVecs() const {return _mcList;}


protected:


private:
  std::vector<Psi2SToKpKmPiGamEvtData*> _dataList;
  std::vector<Psi2SToKpKmPiGamEvtData*> _mcList;

  void read4Vecs(EventList& evtList, std::vector<Psi2SToKpKmPiGamEvtData*>& theEvtList);
  

};

#endif 
