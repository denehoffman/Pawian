#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarp/pbarpDataBaseList.hh"

#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


pbarpDataBaseList::pbarpDataBaseList() :
  EvtDataBaseListNew()
{
}


pbarpDataBaseList::~pbarpDataBaseList()
{
}

void pbarpDataBaseList::read4Vecs(EventList& evtList, std::vector<EvtDataNew*>& theEvtList, double& evtWeightSum, int maxEvts){

  return;
}


