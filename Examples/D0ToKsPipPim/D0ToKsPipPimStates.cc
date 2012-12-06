#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/D0ToKsPipPim/D0ToKsPipPimStates.hh"
#include "ErrLogger/ErrLogger.hh"

D0ToKsPipPimStates::D0ToKsPipPimStates() : 
  _D0JPC(new jpcRes(0, -1, 0)),
  _f0JPC(new jpcRes(0, 1, 1)),
  _piJPC(new jpcRes(0, -1, 1)),
  _kJPC(new jpcRes(0, -1, 0)),
  _k0JPC(new jpcRes(0, 1, 0))
{
  //X decays to Phi Phi
  fillJPCLS(_D0JPC, _f0JPC, _kJPC, _JPCLS_D0ToF0K);
  fillJPCLS(_D0JPC, _k0JPC, _kJPC,_JPCLS_D0ToK0K);

}

D0ToKsPipPimStates::~D0ToKsPipPimStates()
{
}

void D0ToKsPipPimStates::print(std::ostream& os) const
{
  os << "*** D0 -> f0 K:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_D0ToF0K );

  os << "*** D0 -> K0 K:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_D0ToK0K );  
}

void D0ToKsPipPimStates::fillJPCLS(boost::shared_ptr<jpcRes> motherRes, boost::shared_ptr<jpcRes> daughterRes1, boost::shared_ptr<jpcRes> daughterRes2, std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSVec)
{
  // first: check C-parity
  if ( motherRes->C != daughterRes1->C*daughterRes2->C){
    Warning << "C-Parity not valid for the reaction: JPC= " 
            << motherRes->J << " " << motherRes->P << " " << motherRes->C
            << " --> "
            << " JPC= " << daughterRes1->J << " " << daughterRes1->P << " " << daughterRes1->C
            << " and "
            << " JPC= " << daughterRes2->J << " " << daughterRes2->P << " " << daughterRes2->C
            ;  // << endmsg; 
  }
  vector<LS> LSs=GetValidLS(motherRes->J, motherRes->P, daughterRes1->J, daughterRes1->P, daughterRes2->J, daughterRes2->P);

  int num_LS = (int) LSs.size();

  for(int ls = 0; ls < num_LS; ls++){
    Spin L= LSs[ls].L; 
    Spin S= LSs[ls].S;
    int LplusS(L+S);
    if ( LplusS%2 ==0 ){
      boost::shared_ptr<const JPCLS> tmpJPCLS(new JPCLS(motherRes, L, S));
      theJPCLSVec.push_back(tmpJPCLS);
    }
  }
}

void D0ToKsPipPimStates::printDecayJPCLS(std::ostream& os,std::vector< boost::shared_ptr<const JPCLS > > theJPCLS) const{
  
  std::vector< boost::shared_ptr<const JPCLS > >::const_iterator itJPCLS;
  for ( itJPCLS=theJPCLS.begin(); itJPCLS!=theJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }
}
