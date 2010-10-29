#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

Psi2SToKpKmPiGamStates::Psi2SToKpKmPiGamStates() :
  _psiJPC(new jpcRes(1, -1, -1)),
  _chic1JPC(new jpcRes(1, 1, 1)),
  _Kst0JPC(new jpcRes(0, 1, 0)),
  _Kst1JPC(new jpcRes(1, -1, 0)),
  _Kst2JPC(new jpcRes(2, 1, 0)),
  _a0JPC(new jpcRes(0, 1, 1)),
  _a2JPC(new jpcRes(2, 1, 1)),
  _pi0JPC(new jpcRes(0, -1, 1)),
  _kJPC(new jpcRes(0, -1, 0)),
  _gamJPC(new jpcRes(1, -1, -1))
{
  //LS combinations for Psi(2S) decay to Chi_c1 gamma
  fillJPCLS(_psiJPC, _chic1JPC, _gamJPC, _PsiJPCLS);

  //LS combinations for the Chi_c1 decay to K*0-+ K+-
  fillJPCLS(_chic1JPC, _Kst0JPC, _kJPC, _ChiToKst0KJPCLS);

  //LS combinations for the Chi_c1 decay to K*1-+ K+-
  fillJPCLS(_chic1JPC, _Kst1JPC, _kJPC, _ChiToKst1KJPCLS);

  //LS combinations for the Chi_c1 decay to K*2-+ K+-
  fillJPCLS(_chic1JPC, _Kst2JPC, _kJPC, _ChiToKst2KJPCLS);

  //LS combinations for the Chi_c1 decay to a0(980) Pi0
  fillJPCLS(_chic1JPC, _a0JPC, _pi0JPC, _ChiToa0PiJPCLS);

  //LS combinations for the Chi_c1 decay to a0(980) Pi0
  fillJPCLS(_chic1JPC, _a2JPC, _pi0JPC, _ChiToa2PiJPCLS);

  //LS combinations for the K*+- decay to K+- Pi0
  fillJPCLS(_Kst1JPC, _kJPC, _pi0JPC, _Kst1JPCLS);

  //LS combinations for the K*+- decay to K+- Pi0
  fillJPCLS(_Kst2JPC, _kJPC, _pi0JPC, _Kst2JPCLS);

  //LS combinations for the a0(980) decay to K K
  fillJPCLS(_a0JPC, _kJPC, _kJPC, _a0JPCLS);

  //LS combinations for the a0(980) decay to K K
  fillJPCLS(_a2JPC, _kJPC, _kJPC, _a2JPCLS);
}

Psi2SToKpKmPiGamStates::~Psi2SToKpKmPiGamStates()
{
}

void Psi2SToKpKmPiGamStates::fillJPCLS(boost::shared_ptr<jpcRes> motherRes, boost::shared_ptr<jpcRes> daughterRes1, boost::shared_ptr<jpcRes> daughterRes2, std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSVec)
{
  // first: check C-parity
  if ( motherRes->C != daughterRes1->C*daughterRes2->C){
    Warning << "C-Parity not valid for the reaction: JPC= " 
	    << motherRes->J << " " << motherRes->P << " " << motherRes->C
	    << " --> "
	    << " JPC= " << daughterRes1->J << " " << daughterRes1->P << " " << daughterRes1->C
	    << " and "
	    << " JPC= " << daughterRes2->J << " " << daughterRes2->P << " " << daughterRes2->C
	    << endmsg; 
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


void Psi2SToKpKmPiGamStates::print(std::ostream& os) const
{
  os << "*** Psi(2S):  LS combinations for the decay to Psi(2S) gamma*** "<< std::endl;
  
  std::vector< boost::shared_ptr<const JPCLS > >::const_iterator itJPCLS;
  
  for ( itJPCLS=_PsiJPCLS.begin(); itJPCLS!=_PsiJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "*** Chi_c1:  LS combinations for the decay to K*0+- K-+ *** "<< std::endl;
  for ( itJPCLS=_ChiToKst0KJPCLS.begin(); itJPCLS!=_ChiToKst0KJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "*** Chi_c1:  LS combinations for the decay to K*1+- K-+ *** "<< std::endl;
  for ( itJPCLS=_ChiToKst1KJPCLS.begin(); itJPCLS!=_ChiToKst1KJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "*** Chi_c1:  LS combinations for the decay to K*2+- K-+ *** "<< std::endl;
  for ( itJPCLS=_ChiToKst2KJPCLS.begin(); itJPCLS!=_ChiToKst2KJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "*** Chi_c1:  LS combinations for the decay to a0(980) pi0 *** "<< std::endl;
  for ( itJPCLS=_ChiToa0PiJPCLS.begin(); itJPCLS!=_ChiToa0PiJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "*** Chi_c1:  LS combinations for the decay to a2 pi0 *** "<< std::endl;
  for ( itJPCLS=_ChiToa2PiJPCLS.begin(); itJPCLS!=_ChiToa2PiJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** K*1+-:   LS combinations for the decay to K*1+- -> K+- pi0 *** "<< std::endl;
  for ( itJPCLS=_Kst1JPCLS.begin(); itJPCLS!=_Kst1JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** K*2+-:   LS combinations for the decay to K*2+- -> K+- pi0 *** "<< std::endl;
  for ( itJPCLS=_Kst2JPCLS.begin(); itJPCLS!=_Kst2JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** a0(980):  LS combinations for the decay to K K *** "<< std::endl;
  for ( itJPCLS=_a0JPCLS.begin(); itJPCLS!=_a0JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** a2:  LS combinations for the decay to K K *** "<< std::endl;
  for ( itJPCLS=_a2JPCLS.begin(); itJPCLS!=_a2JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** K+- *** "<< std::endl;
  _kJPC->print(os);

  os << "\n *** pi0 *** "<< std::endl;
  _pi0JPC->print(os);

  os << "\n *** gamma *** "<< std::endl;
  _gamJPC->print(os);
  os << "\n"<< std::endl;
}
