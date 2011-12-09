#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKStates.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamKsKlKKStates::JpsiGamKsKlKKStates() :
  _psiJPC(new jpcRes(1, -1, -1)),
  _gammaJPC(new jpcRes(1, -1, -1)),
  _kJPC(new jpcRes(0, -1, 0)),
  _phiJPC(new jpcRes(1, -1, -1)),
  _etaJPC(new jpcRes(0, -1, 1)),
  _f0JPC(new jpcRes(0, 1, 1)),
  _f2JPC(new jpcRes(2, 1, 1))
{
  //LS combinations for Psi decays to X gamma
  fillJPCLS(_psiJPC, _etaJPC, _gammaJPC, _JPCLS_PsiToEtaGamma);
  fillJPCLS(_psiJPC, _f0JPC, _gammaJPC, _JPCLS_PsiToF0Gamma);
  fillJPCLS(_psiJPC, _f2JPC, _gammaJPC, _JPCLS_PsiToF2Gamma);
  
  //X decays to phi phi
  fillJPCLS(_etaJPC, _phiJPC, _phiJPC, _JPCLS_EtaToPhiPhi);
  fillJPCLS(_f0JPC, _phiJPC, _phiJPC, _JPCLS_F0ToPhiPhi);
  fillJPCLS(_f2JPC, _phiJPC, _phiJPC, _JPCLS_F2ToPhiPhi);

  //phi to kk
  fillJPCLS(_phiJPC, _kJPC, _kJPC, _JPCLS_PhiToKK);
  
}

JpsiGamKsKlKKStates::~JpsiGamKsKlKKStates()
{
}

void JpsiGamKsKlKKStates::fillJPCLS(boost::shared_ptr<jpcRes> motherRes, boost::shared_ptr<jpcRes> daughterRes1, boost::shared_ptr<jpcRes> daughterRes2, std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSVec)
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
    boost::shared_ptr<const JPCLS> tmpJPCLS(new JPCLS(motherRes, L, S));
    theJPCLSVec.push_back(tmpJPCLS);
  }
} 


void JpsiGamKsKlKKStates::print(std::ostream& os) const
{
  os << "*** Psi -> Eta gamma:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_PsiToEtaGamma );
  
  os << "*** Psi -> f0 gamma:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_PsiToF0Gamma );
  
  os << "*** Psi -> f2 gamma:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_PsiToF2Gamma );
  
  os << "*** eta -> phi phi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_EtaToPhiPhi );
  
  os << "*** f0 -> phi phi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_F0ToPhiPhi );
  
  os << "*** f2 -> phi phi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_F2ToPhiPhi );
  
  os << "*** phi -> KK:   LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_PhiToKK );
}


void JpsiGamKsKlKKStates::printDecayJPCLS(std::ostream& os,std::vector< boost::shared_ptr<const JPCLS > > theJPCLS) const{
  
  std::vector< boost::shared_ptr<const JPCLS > >::const_iterator itJPCLS;
  for ( itJPCLS=theJPCLS.begin(); itJPCLS!=theJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }
}
