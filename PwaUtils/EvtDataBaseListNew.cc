#include <getopt.h>


#include "PwaUtils/EvtDataBaseListNew.hh"
#include "Event/EventList.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


EvtDataBaseListNew::EvtDataBaseListNew() :
  _noOfWeightedDataEvts(0.),
  _noOfWeightedMcEvts(0.),
  _mcToDataRatio(1000),
  _alreadyRead(false)
{
}


EvtDataBaseListNew::~EvtDataBaseListNew()
{
}

void EvtDataBaseListNew::read(EventList& evtListData, EventList& evtListMc){
  if(_alreadyRead){
   Alert << "4 vectors already read " << endmsg;  // << endmsg;
    exit(1);
  }
  read4Vecs(evtListData, _evtDataList, _noOfWeightedDataEvts, evtListData.size());

  int maxMcEvts=evtListMc.size();
  if (maxMcEvts > _mcToDataRatio*evtListData.size() ) maxMcEvts=_mcToDataRatio*evtListData.size();
  read4Vecs(evtListMc, _mcDataList, _noOfWeightedMcEvts, maxMcEvts);
  _alreadyRead=true;
}

