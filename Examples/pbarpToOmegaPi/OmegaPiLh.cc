#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPi/OmegaPiLh.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "PwaUtils/pbarpStates.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"


OmegaPiLh::OmegaPiLh(boost::shared_ptr<const OmegaPiEventList> theEvtList, int jmax) :
  _omegaPiEventListPtr(theEvtList),
  _jmax(jmax)
{
  init();
}

OmegaPiLh::OmegaPiLh(boost::shared_ptr<OmegaPiLh> theOmegaPiLhPtr):
  _omegaPiEventListPtr(theOmegaPiLhPtr->getEventList()),
  _jmax(theOmegaPiLhPtr->jmax())
{
  init();
}

OmegaPiLh::~OmegaPiLh()
{
}


void OmegaPiLh::init(){
  boost::shared_ptr<pbarpStates> thepbarpStatesPtr(new pbarpStates(_jmax));
  boost::shared_ptr<pbarpToOmegaPi0States> theOmegaPi0StatesPtr(new pbarpToOmegaPi0States(thepbarpStatesPtr));
   _JPCLSomegaProd=theOmegaPi0StatesPtr->jpclsStates();
   _pbarpJPC=theOmegaPi0StatesPtr->jpcStates();
   _pbarpSinglett=thepbarpStatesPtr->extractSingletStates(_pbarpJPC);
}

void OmegaPiLh::print(std::ostream& os) const{
  os << "OmegaPiLh::print\n"
     << "The contributed JPC initial states are:" << std::endl;
  
  std::vector< boost::shared_ptr<const jpcRes > >::const_iterator itJPC;
  for ( itJPC=_pbarpJPC.begin(); itJPC!=_pbarpJPC.end(); ++itJPC){
    (*itJPC)->print(os);
    os << std::endl;
  }

  os << "The initial JPC singlett states are:" << std::endl;
  for ( itJPC=_pbarpSinglett.begin(); itJPC!=_pbarpSinglett.end(); ++itJPC){
    (*itJPC)->print(os);
    os << std::endl;
  }

  os << "The contributed JPC LS combinations for the omega production are:" << std::endl;
  std::vector< boost::shared_ptr<const JPCLS > >::const_iterator itJPCLS;
  for ( itJPCLS=_JPCLSomegaProd.begin(); itJPCLS!=_JPCLSomegaProd.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }
  
}
