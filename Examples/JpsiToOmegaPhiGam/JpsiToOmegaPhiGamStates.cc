#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/JpsiToOmegaPhiGam/JpsiToOmegaPhiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiToOmegaPhiGamStates::JpsiToOmegaPhiGamStates() : 
  PsiToXGamStates(),  
  _kJPC(new jpcRes(0, -1, 0)),
  _phiJPC(new jpcRes(1, -1, -1)),
  _etaJPC(new jpcRes(0, -1, 1)),
  _f0JPC(new jpcRes(0, 1, 1)),
  _f1JPC(new jpcRes(1, 1, 1)),
  _f2JPC(new jpcRes(2, 1, 1)),
  _eta2JPC(new jpcRes(2,-1,1) )
{
  //X decays to Phi Phi
  fillJPCLS(_etaJPC, _phiJPC, _phiJPC, _JPCLS_EtaToPhiPhi);
  fillJPCLS(_f0JPC, _phiJPC, _phiJPC, _JPCLS_F0ToPhiPhi);
  fillJPCLS(_f1JPC, _phiJPC, _phiJPC, _JPCLS_F1ToPhiPhi);
  fillJPCLS(_f2JPC, _phiJPC, _phiJPC, _JPCLS_F2ToPhiPhi);
  fillJPCLS(_eta2JPC, _phiJPC, _phiJPC, _JPCLS_Eta2ToPhiPhi);
  
  //phi to kk
  fillJPCLS(_phiJPC, _kJPC, _kJPC, _JPCLS_PhiToKK);

}

JpsiToOmegaPhiGamStates::~JpsiToOmegaPhiGamStates()
{
}

void JpsiToOmegaPhiGamStates::print(std::ostream& os) const
{
  PsiToXGamStates::print(os);

  os << "*** eta -> phi phi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_EtaToPhiPhi );
  
  os << "*** f0 -> phi phi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_F0ToPhiPhi );
  
  os << "*** f1 -> phi phi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_F1ToPhiPhi );

  os << "*** f2 -> phi phi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_F2ToPhiPhi );
  
  os << "*** eta2 -> phi phi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_Eta2ToPhiPhi );
  
  os << "*** phi -> KK:   LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_PhiToKK );
}

void JpsiToOmegaPhiGamStates::fillJPCLS(boost::shared_ptr<jpcRes> motherRes, boost::shared_ptr<jpcRes> daughterRes1, boost::shared_ptr<jpcRes> daughterRes2, std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSVec)
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
