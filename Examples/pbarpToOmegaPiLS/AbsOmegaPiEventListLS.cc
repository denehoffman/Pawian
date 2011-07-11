#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"
#include "Event/EventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "Event/Event.hh"
#include "ErrLogger/ErrLogger.hh"


AbsOmegaPiEventListLS::AbsOmegaPiEventListLS(EventList& evtListData, EventList& evtListMc, unsigned jmax, unsigned pbarmom):
  _jmax(jmax),
  _pbarmom(pbarmom)
{
  if (_jmax<0) { Alert << "_lmax < 0 is not allowed!!!" << endmsg; exit(1); }
}




AbsOmegaPiEventListLS::~AbsOmegaPiEventListLS()
{
}



