#ifndef _Psi2STo2K2PiGamEvtList_H
#define _Psi2STo2K2PiGamEvtList_H

#include <iostream>
#include <vector>

#include <cassert>
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"

using Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData;

class EventList;
class Event;

class Psi2STo2K2PiGamEvtList {

public:

  // create/copy/destroy:

  ///Constructor  
  Psi2STo2K2PiGamEvtList(); 
  Psi2STo2K2PiGamEvtList(EventList& evtListData, EventList& evtListMc);



  /** Destructor */
  virtual ~Psi2STo2K2PiGamEvtList();

  Psi2STo2K2PiGamEvtData* fillEvtData(Event* anEvent, int evtNo=0);
  // Getters (not working for Psi2STo2K2PiGamEvtList() constructor:
  const std::vector<Psi2STo2K2PiGamEvtData*> getDataVecs() const {return _dataList;}
  const std::vector<Psi2STo2K2PiGamEvtData*> getMcVecs() const {return _mcList;}


protected:


private:
  std::vector<Psi2STo2K2PiGamEvtData*> _dataList;
  std::vector<Psi2STo2K2PiGamEvtData*> _mcList;

  void read4Vecs(EventList& evtList, std::vector<Psi2STo2K2PiGamEvtData*>& theEvtList);
  

};

#endif 
