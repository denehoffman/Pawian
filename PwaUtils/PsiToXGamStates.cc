#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "PwaUtils/PsiToXGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

PsiToXGamStates::PsiToXGamStates() :
  _psiJPC(new jpcRes(1, -1, -1)),
  _gammaJPC(new jpcRes(1, -1, -1)),
  _etaJPC(new jpcRes(0, -1, 1)),
  _f0JPC(new jpcRes(0, 1, 1)),
  _f1JPC(new jpcRes(1, 1, 1)),
  _f2JPC(new jpcRes(2, 1, 1)),
  _eta2JPC(new jpcRes(2,-1,1) )
{
  //LS combinations for Psi decays to X gamma
  fillJPCLS(_psiJPC, _etaJPC , _gammaJPC, _JPCLS_PsiToEtaGamma );
  fillJPCLS(_psiJPC, _f0JPC  , _gammaJPC, _JPCLS_PsiToF0Gamma  );
  fillJPCLS(_psiJPC, _f1JPC  , _gammaJPC, _JPCLS_PsiToF1Gamma  );
  fillJPCLS(_psiJPC, _f2JPC  , _gammaJPC, _JPCLS_PsiToF2Gamma  );
  fillJPCLS(_psiJPC, _eta2JPC, _gammaJPC, _JPCLS_PsiToEta2Gamma);

  fillJPClamlam(_JPCLS_PsiToEtaGamma,  _etaJPC,  _JPCLamLam_PsiToEtaGamma);
  fillJPClamlam(_JPCLS_PsiToF0Gamma,   _f0JPC,   _JPCLamLam_PsiToF0Gamma);
  fillJPClamlam(_JPCLS_PsiToF2Gamma,   _f1JPC,   _JPCLamLam_PsiToF1Gamma);
  fillJPClamlam(_JPCLS_PsiToF2Gamma,   _f2JPC,   _JPCLamLam_PsiToF2Gamma);
  fillJPClamlam(_JPCLS_PsiToEta2Gamma, _eta2JPC, _JPCLamLam_PsiToEta2Gamma);  
}

PsiToXGamStates::~PsiToXGamStates()
{
}

void PsiToXGamStates::fillJPCLS(boost::shared_ptr<jpcRes> motherRes, boost::shared_ptr<jpcRes> daughterRes1, boost::shared_ptr<jpcRes> daughterRes2, std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSVec)
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

void PsiToXGamStates::fillJPClamlam(std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSvec, boost::shared_ptr<jpcRes> xRes, std::vector< boost::shared_ptr<const JPClamlam> >& toFill){
  //retrieve Smax
  Spin Smax(0);
  
  std::vector< boost::shared_ptr<const JPCLS> >::iterator it;

  for (it=theJPCLSvec.begin(); it!=theJPCLSvec.end(); ++it){
    Spin currentS=(*it)->S;
    if ( Smax < currentS) Smax=currentS;
  }
  Info << "Smax= " << Smax << endmsg;

  Spin  lamgam(1);
  for (Spin lamx = 0; lamx <= xRes->J; ++lamx){
      if (fabs(lamx-lamgam)>_psiJPC->J || fabs(lamx-lamgam)>Smax) continue;
      double parityFactor = _psiJPC->P * xRes->P * _gammaJPC->P * pow(-1., xRes->J+_gammaJPC->J-_psiJPC->J );
      boost::shared_ptr<const JPClamlam> tmpJPClamlam(new JPClamlam(_psiJPC, lamx, lamgam, parityFactor));
      toFill.push_back(tmpJPClamlam);
  }
  
}

void PsiToXGamStates::print(std::ostream& os) const
{
  os << "*** Psi -> Eta gamma:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_PsiToEtaGamma );

  os << "*** Psi -> Eta gamma:  lambda lambda combinations for the decay *** "<< std::endl;
  printDecayJPClamlam(os, _JPCLamLam_PsiToEtaGamma );

  os << "*** Psi -> f0 gamma:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_PsiToF0Gamma );

  os << "*** Psi -> f0 gamma:  lambda lambda combinations for the decay *** "<< std::endl;
  printDecayJPClamlam(os, _JPCLamLam_PsiToF0Gamma);  

  os << "*** Psi -> f1 gamma:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_PsiToF1Gamma );

  os << "*** Psi -> f1 gamma:  lambda lambda combinations for the decay *** "<< std::endl;
  printDecayJPClamlam(os, _JPCLamLam_PsiToF1Gamma); 

  os << "*** Psi -> f2 gamma:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_PsiToF2Gamma );

  os << "*** Psi -> f2 gamma:  lambda lambda combinations for the decay *** "<< std::endl;
  printDecayJPClamlam(os, _JPCLamLam_PsiToF2Gamma);
  
  os << "*** Psi -> eta2 gamma:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_PsiToEta2Gamma );

  os << "*** Psi -> eta2 gamma:  lambda lambda  combinations for the decay *** "<< std::endl;
  printDecayJPClamlam(os, _JPCLamLam_PsiToEta2Gamma );
}


void PsiToXGamStates::printDecayJPCLS(std::ostream& os,std::vector< boost::shared_ptr<const JPCLS > > theJPCLS) const{
  
  std::vector< boost::shared_ptr<const JPCLS > >::const_iterator itJPCLS;
  for ( itJPCLS=theJPCLS.begin(); itJPCLS!=theJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << "\n" << std::endl;
  }
}
 

void PsiToXGamStates::printDecayJPClamlam(std::ostream& os,std::vector< boost::shared_ptr<const JPClamlam > > theJPClamlam) const{

  std::vector< boost::shared_ptr<const JPClamlam > >::const_iterator it;
  for ( it=theJPClamlam.begin(); it!=theJPClamlam.end(); ++it){
    (*it)->print(os);
    os << "\n" << std::endl;
  }
}
