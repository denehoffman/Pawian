#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPi/OmegaPiLh.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "PwaUtils/pbarpStates.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"


OmegaPiLh::OmegaPiLh(boost::shared_ptr<const OmegaPiEventList> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0States> theStates) :
  _omegaPiEventListPtr(theEvtList),
  _omegaPi0StatesPtr(theStates)
{
}

OmegaPiLh::OmegaPiLh(boost::shared_ptr<OmegaPiLh> theOmegaPiLhPtr):
  _omegaPiEventListPtr(theOmegaPiLhPtr->getEventList()),
  _omegaPi0StatesPtr(theOmegaPiLhPtr->omegaPi0States())
{
}

OmegaPiLh::~OmegaPiLh()
{
}


void OmegaPiLh::print(std::ostream& os) const{
  os << "OmegaPiLh::print\n";
}
