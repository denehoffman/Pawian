#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/Tutorial/DfuncClebschG/JPCdecays.hh"
#include "PwaUtils/DataUtils.hh"



#include "ErrLogger/ErrLogger.hh"

JPCdecays::JPCdecays(boost::shared_ptr<const jpcRes> motherJPC, boost::shared_ptr<const jpcRes> daughter1JPC, boost::shared_ptr<const jpcRes> daughter2JPC, bool sameDaughters) :
  _motherJPC(motherJPC)
  ,_daughter1JPC(daughter1JPC)
  ,_daughter2JPC(daughter2JPC)
  ,_sameDaughters(sameDaughters)
{
  _allLSs=GetValidLS(_motherJPC->J, _motherJPC->P, _daughter1JPC->J,_daughter1JPC->P, _daughter2JPC->J, _daughter2JPC->P);


  //check C parity
  bool allowedC=true;
  if (!sameDaughters){
    if (_motherJPC->C != (daughter1JPC->C * daughter2JPC->C)) allowedC=false;
  }

  std::vector<LS>::const_iterator itLS;

  for (itLS=_allLSs.begin(); itLS!=_allLSs.end(); ++itLS){
    if (sameDaughters){
      bool bosons=true;
      if (int(2*daughter1JPC->J)%2 !=0) bosons=false;
      
      Spin theL=itLS->L;
      Spin theS=itLS->S;
      int cparity=0;

      if (bosons) cparity=int(pow(-1,theL));
      else cparity=int(pow(-1,(theL+theS))); //daughters are fermions	

      if(_motherJPC->C != cparity ) allowedC=false;
    }

    if (allowedC) _allowedCParityLSs.push_back(*itLS);
    else _forbiddenCParityLSs.push_back(*itLS);
  }
  
}

JPCdecays::~JPCdecays()
{

}

void JPCdecays::print(std::ostream& os) const{

  std::cout << "the following LS combinations are allowed for the decay: ";
  _motherJPC->print(os);
  std::cout << " -> ";
  _daughter1JPC->print(os);
  std::cout << " + ";
  _daughter2JPC->print(os);

  std::cout << "\n __________LS combinations with allowed C parity______________" << std::endl;

  std::vector<LS>::const_iterator itLS;

  for (itLS=_allowedCParityLSs.begin(); itLS!=_allowedCParityLSs.end(); ++itLS){
    std::cout << "L: " << itLS->L << "\tS: " << itLS->S << std::endl;
  }

 std::cout << "\n __________LS combinations with forbidden C parity______________" << std::endl;

  for (itLS=_forbiddenCParityLSs.begin(); itLS!=_forbiddenCParityLSs.end(); ++itLS){
    std::cout << "L: " << itLS->L << "\tS: " << itLS->S << std::endl;
  }
}
