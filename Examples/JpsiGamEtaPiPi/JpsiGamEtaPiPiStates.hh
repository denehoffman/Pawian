#ifndef _JpsiGamEtaPiPiStates_H
#define _JpsiGamEtaPiPiStates_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>


#include <boost/shared_ptr.hpp>

// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"


class JpsiGamEtaPiPiStates {

public:

  ///Constructor 
  JpsiGamEtaPiPiStates();


  /** Destructor */
  virtual ~JpsiGamEtaPiPiStates();

  std::vector< boost::shared_ptr<const JPCLS> > PsiToEtaGammaStates() const {return _JPCLS_PsiToEtaGamma;}
  std::vector< boost::shared_ptr<const JPCLS> > PsiToF0GammaStates() const {return _JPCLS_PsiToF0Gamma;}
  std::vector< boost::shared_ptr<const JPCLS> > PsiToF1GammaStates() const {return _JPCLS_PsiToF1Gamma;}
  std::vector< boost::shared_ptr<const JPCLS> > PsiToF2GammaStates() const {return _JPCLS_PsiToF2Gamma;}
  /*
  std::vector< boost::shared_ptr<const JPCLS> > PsiToEta2GammaStates() const {return _JPCLS_PsiToEta2Gamma;}
  std::vector< boost::shared_ptr<const JPCLS> > EtaToPhiPhiStates() const {return _JPCLS_EtaToPhiPhi;}
  std::vector< boost::shared_ptr<const JPCLS> > F0ToPhiPhiStates() const {return _JPCLS_F0ToPhiPhi;}
  std::vector< boost::shared_ptr<const JPCLS> > F1ToPhiPhiStates() const {return _JPCLS_F1ToPhiPhi;}
  std::vector< boost::shared_ptr<const JPCLS> > F2ToPhiPhiStates() const {return _JPCLS_F2ToPhiPhi;}
  std::vector< boost::shared_ptr<const JPCLS> > Eta2ToPhiPhiStates() const {return _JPCLS_Eta2ToPhiPhi;}
  std::vector< boost::shared_ptr<const JPCLS> > PhiToKKStates() const {return _JPCLS_PhiToKK;}
  */
  
  void print(std::ostream& os) const; 

protected:


private:

  boost::shared_ptr<jpcRes> _psiJPC;
  boost::shared_ptr<jpcRes> _gammaJPC;
  boost::shared_ptr<jpcRes> _etaJPC;
  boost::shared_ptr<jpcRes> _f0JPC;
  boost::shared_ptr<jpcRes> _f1JPC;
  boost::shared_ptr<jpcRes> _f2JPC;
  /*
  boost::shared_ptr<jpcRes> _kJPC;
  boost::shared_ptr<jpcRes> _phiJPC;
  boost::shared_ptr<jpcRes> _eta2JPC;
  */
  
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_PsiToEtaGamma;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_PsiToF0Gamma;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_PsiToF1Gamma;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_PsiToF2Gamma;
  /*
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_PsiToEta2Gamma;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_EtaToPhiPhi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_F0ToPhiPhi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_F1ToPhiPhi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_F2ToPhiPhi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_Eta2ToPhiPhi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_PhiToKK;
  */
  
  void fillJPCLS(boost::shared_ptr<jpcRes>, boost::shared_ptr<jpcRes>, boost::shared_ptr<jpcRes>, std::vector< boost::shared_ptr<const JPCLS> >& );
  void printDecayJPCLS(std::ostream& os,std::vector< boost::shared_ptr<const JPCLS > > theJPCLS) const;
  
};

#endif
