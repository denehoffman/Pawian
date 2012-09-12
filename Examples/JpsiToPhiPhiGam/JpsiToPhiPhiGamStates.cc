#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/JpsiToPhiPhiGam/JpsiToPhiPhiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiToPhiPhiGamStates::JpsiToPhiPhiGamStates() : 
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

JpsiToPhiPhiGamStates::~JpsiToPhiPhiGamStates()
{
}

void JpsiToPhiPhiGamStates::print(std::ostream& os) const
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


