#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

Psi2STo2K2PiGamStates::Psi2STo2K2PiGamStates() :
  _psiJPC(new jpcRes(1, -1, -1)),
  _chic0JPC(new jpcRes(0, 1, 1)),
  _Kst0JPC(new jpcRes(0, 1, 0)),
  _Kst1JPC(new jpcRes(1, -1, 0)),
  _Kst2JPC(new jpcRes(2, 1, 0)),
  _K1400JPC(new jpcRes(1, 1, 0)),
  _f0JPC(new jpcRes(0, 1, 1)),
  _f2JPC(new jpcRes(2, 1, 1)),
  _pi0JPC(new jpcRes(0, -1, 1)),
  _kJPC(new jpcRes(0, -1, 0)),
  _gamJPC(new jpcRes(1, -1, -1)),
  _pi2JPC(new jpcRes(2, -1, 1)),
  _K2mJPC(new jpcRes(2, -1, 0))
{
  //LS combinations for Psi(2S) decay to Chi_c0 gamma
  fillJPCLS(_psiJPC, _chic0JPC, _gamJPC, _PsiJPCLS);

  //LS combinations for the Chi_c0 decay to K1(1400) K
  fillJPCLS(_chic0JPC, _K1400JPC, _kJPC, _ChiToK1400ToK892piJPCLS);

  //LS combinations for the Chi_c0 decay to K*(892) K*(892)
  fillJPCLS(_chic0JPC, _Kst1JPC, _Kst1JPC,_ChiTo2K892JPCLS);

  //LS combinations for the Chi_c0 decay to K*2(1430) K*2(1430)
  fillJPCLS(_chic0JPC, _Kst2JPC, _Kst2JPC, _ChiTo2K_2_1430JPCLS);

  //LS combinations for the Chi_c0 decay to K*0(1430) K*0(1430)
  fillJPCLS(_chic0JPC, _Kst0JPC, _Kst0JPC, _ChiTo2K_0_JPCLS);

  //LS combinations for the Chi_c0 decay to K*0(1430) K*2(1430)
  fillJPCLS(_chic0JPC, _Kst0JPC, _Kst2JPC, _ChiToK0K2_JPCLS);

  //LS combinations for the Chi_c0 decay to f0 f0
  fillJPCLS(_chic0JPC, _f0JPC, _f0JPC,_ChiTof0f0JPCLS);

  //LS combinations for the Chi_c0 decay to f0 f2
  fillJPCLS(_chic0JPC, _f0JPC, _f2JPC,_ChiTof0f2JPCLS);

  //LS combinations for the Chi_c0 decay to K0 K0
  fillJPCLS(_chic0JPC, _kJPC, _kJPC, _ChiToK0K0JPCLS);

  //LS combinations for the Chi_c0 decay to Pi2 Pi
  fillJPCLS(_chic0JPC, _pi0JPC, _pi2JPC, _ChiToPi_2PiJPCLS);

  //LS combinations for the K*+- decay to K+- Pi0
  fillJPCLS(_Kst1JPC, _kJPC, _pi0JPC, _Kst1JPCLS);

  //LS combinations for the K*2+- decay to K+- Pi0
  fillJPCLS(_Kst2JPC, _kJPC, _pi0JPC, _Kst2JPCLS);

  //LS combinations for the K*+- decay to K+- Pi0
  fillJPCLS(_K1400JPC, _Kst1JPC, _pi0JPC, _K1400ToKst1PiJPCLS);

  //LS combinations for the K1- decay to K0 Pi0
  fillJPCLS(_K1400JPC, _Kst0JPC, _pi0JPC, _K1ToK0PiJPCLS);

  //LS combinations for the f2 decay to K+ K-
  fillJPCLS(_f2JPC, _kJPC, _kJPC, _f2JPCLS); 

 //LS combinations for the pi_2 decay to f2 Pi
  fillJPCLS(_pi2JPC, _pi0JPC, _f2JPC, _Pi_2Tof_2PiJPCLS);

 //LS combinations for the pi_2 decay to K*0 K 
  fillJPCLS(_pi2JPC, _kJPC, _Kst0JPC,_Pi_2ToKst0KJPCLS);

 //LS combinations for the pi_2 decay to K*1 K
  fillJPCLS(_pi2JPC, _kJPC, _Kst1JPC,_Pi_2ToKst1KJPCLS);

 //LS combinations for the pi_2 decay to K*2 K 
  fillJPCLS(_pi2JPC, _kJPC, _Kst2JPC,_Pi_2ToKst2KJPCLS);

 //LS combinations for the chi_c0 to pi0 pi0
  fillJPCLS(_chic0JPC, _pi0JPC, _pi0JPC, _ChiToPi0Pi0JPCLS);

 //LS combinations for the pi_2 decay to f0 Pi
  fillJPCLS(_pi2JPC, _pi0JPC, _f0JPC, _Pi_2Tof0PiJPCLS);

 //LS combinations for the decay chi_c0 to K2- K0-
  fillJPCLS(_chic0JPC, _kJPC, _K2mJPC, _ChiToK2mK0mJPCLS);

 //LS combinations for the decay K2- pi0
  fillJPCLS(_K2mJPC, _pi0JPC, _Kst2JPC,_K2mToK2pPiJPCLS);

 //LS combinations for the K1+ decay to f0 K
  fillJPCLS(_K1400JPC, _f0JPC, _kJPC, _K1pTof0KJPCLS);

  fillJPCLS(_pi0JPC, _Kst1JPC, _pi0JPC, _Pi0pToKstKJPCLS);

}

Psi2STo2K2PiGamStates::~Psi2STo2K2PiGamStates()
{
}

void Psi2STo2K2PiGamStates::fillJPCLS(boost::shared_ptr<jpcRes> motherRes, boost::shared_ptr<jpcRes> daughterRes1, boost::shared_ptr<jpcRes> daughterRes2, std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSVec)
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


void Psi2STo2K2PiGamStates::print(std::ostream& os) const
{
  os << "*** Psi(2S):  LS combinations for the decay to Psi(2S) gamma*** "<< std::endl;
  
  std::vector< boost::shared_ptr<const JPCLS > >::const_iterator itJPCLS;
  
  for ( itJPCLS=_PsiJPCLS.begin(); itJPCLS!=_PsiJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "*** Chi_c0:  LS combinations for the decay to K1(1400)*+- K-+ *** "<< std::endl;
  for ( itJPCLS=_ChiToK1400ToK892piJPCLS.begin(); itJPCLS!=_ChiToK1400ToK892piJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "*** Chi_c0:  LS combinations for the decay to K*(892) K*(892) *** "<< std::endl;
  for ( itJPCLS=_ChiTo2K892JPCLS.begin(); itJPCLS!=_ChiTo2K892JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }


  os << "*** Chi_c0:  LS combinations for the decay to K*2(1430) K*2(1430) *** "<< std::endl;
  for ( itJPCLS=_ChiTo2K_2_1430JPCLS.begin(); itJPCLS!=_ChiTo2K_2_1430JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "*** Chi_c0:  LS combinations for the decay to K*0(1430) K*0(1430) *** "<< std::endl;
  for ( itJPCLS=_ChiTo2K_0_JPCLS.begin(); itJPCLS!=_ChiTo2K_0_JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "*** Chi_c0:  LS combinations for the decay to K*0(1430) K*2(1430) *** "<< std::endl;
  for ( itJPCLS=_ChiToK0K2_JPCLS.begin(); itJPCLS!=_ChiToK0K2_JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "*** Chi_c0:  LS combinations for the decay to f0 f0 *** "<< std::endl;
  for ( itJPCLS=_ChiTof0f0JPCLS.begin(); itJPCLS!=_ChiTof0f0JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "*** Chi_c0:  LS combinations for the decay to f0 f2 *** "<< std::endl;
  for ( itJPCLS=_ChiTof0f2JPCLS.begin(); itJPCLS!=_ChiTof0f2JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "*** Chi_c0:  LS combinations for the decay to pi2 pi0 *** "<< std::endl;
  for ( itJPCLS=_ChiToPi_2PiJPCLS.begin(); itJPCLS!=_ChiToPi_2PiJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** K*2+-:   LS combinations for the decay to K*2+- -> K+- pi0 *** "<< std::endl;
  for ( itJPCLS=_Kst2JPCLS.begin(); itJPCLS!=_Kst2JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** K1(1400):   LS combinations for the decay to K1(1400)+- -> K*1+- pi0 *** "<< std::endl;
  for ( itJPCLS=_K1400ToKst1PiJPCLS.begin(); itJPCLS!=_K1400ToKst1PiJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** K1(1270):   LS combinations for the decay to K1(1270) -> K0(1430) pi0 *** "<< std::endl;
  for ( itJPCLS=_K1ToK0PiJPCLS.begin(); itJPCLS!=_K1ToK0PiJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** f2:  LS combinations for the decay to K K *** "<< std::endl;
  for ( itJPCLS=_f2JPCLS.begin(); itJPCLS!=_f2JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** pi2:  LS combinations for the decay to K0 K0 *** "<< std::endl;
  for ( itJPCLS=_ChiToK0K0JPCLS.begin(); itJPCLS!=_ChiToK0K0JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** pi2:  LS combinations for the pi2 decay to pi0 f2 *** "<< std::endl;
  for ( itJPCLS=_Pi_2Tof_2PiJPCLS.begin(); itJPCLS!=_Pi_2Tof_2PiJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** pi2:  LS combinations for the pi2 decay to K*0 K *** "<< std::endl;
  for ( itJPCLS=_Pi_2ToKst0KJPCLS.begin(); itJPCLS!=_Pi_2ToKst0KJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** pi2:  LS combinations for the pi2 decay to K*1 K *** "<< std::endl;
  for ( itJPCLS=_Pi_2ToKst1KJPCLS.begin(); itJPCLS!=_Pi_2ToKst1KJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** pi2:  LS combinations for the pi2 decay to K*2 K *** "<< std::endl;
  for ( itJPCLS=_Pi_2ToKst2KJPCLS.begin(); itJPCLS!=_Pi_2ToKst2KJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** pi2:  LS combinations for the chi_c0 decay to pi0 pi0 *** "<< std::endl;
  for ( itJPCLS=_ChiToPi0Pi0JPCLS.begin(); itJPCLS!=_ChiToPi0Pi0JPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }


  os << "\n *** pi2:  LS combinations for the pi2 decay to f0 pi0 *** "<< std::endl;
  for ( itJPCLS=_Pi_2Tof0PiJPCLS.begin(); itJPCLS!=_Pi_2Tof0PiJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** pi2:  LS combinations for the chi_c0 decay to K2- K0- *** "<< std::endl;
  for ( itJPCLS=_ChiToK2mK0mJPCLS.begin(); itJPCLS!=_ChiToK2mK0mJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** pi2:  LS combinations for the decay K2- to K*20+ pi0 *** "<< std::endl;
  for ( itJPCLS=_K2mToK2pPiJPCLS.begin(); itJPCLS!=_K2mToK2pPiJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n *** K1+:  LS combinations for the decay K1+ to f0 K *** "<< std::endl;
  for ( itJPCLS=_K1pTof0KJPCLS.begin(); itJPCLS!=_K1pTof0KJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }

  os << "\n ***  LS combinations for the decay pi0 to K*1 K *** "<< std::endl;
  for ( itJPCLS=_Pi0pToKstKJPCLS.begin(); itJPCLS!=_Pi0pToKstKJPCLS.end(); ++itJPCLS){
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
