#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>


#include <boost/shared_ptr.hpp>

// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/PsiToXGamStates.hh"


class JpsiToOmegaPhiGamStates : public PsiToXGamStates {

public:

  ///Constructor 
  JpsiToOmegaPhiGamStates();

  /** Destructor */
  virtual ~JpsiToOmegaPhiGamStates();

  std::vector< boost::shared_ptr<const JPCLS> > EtaToPhiPhiStates() const {return _JPCLS_EtaToPhiPhi;}
  std::vector< boost::shared_ptr<const JPCLS> > F0ToPhiPhiStates() const {return _JPCLS_F0ToPhiPhi;}
  std::vector< boost::shared_ptr<const JPCLS> > F1ToPhiPhiStates() const {return _JPCLS_F1ToPhiPhi;}
  std::vector< boost::shared_ptr<const JPCLS> > F2ToPhiPhiStates() const {return _JPCLS_F2ToPhiPhi;}
  std::vector< boost::shared_ptr<const JPCLS> > Eta2ToPhiPhiStates() const {return _JPCLS_Eta2ToPhiPhi;}
  std::vector< boost::shared_ptr<const JPCLS> > PhiToKKStates() const {return _JPCLS_PhiToKK;}
 
  virtual void print(std::ostream& os) const; 

protected:

private:
  boost::shared_ptr<jpcRes> _kJPC;
  boost::shared_ptr<jpcRes> _phiJPC;
  boost::shared_ptr<jpcRes> _etaJPC;
  boost::shared_ptr<jpcRes> _f0JPC;
  boost::shared_ptr<jpcRes> _f1JPC;
  boost::shared_ptr<jpcRes> _f2JPC;
  boost::shared_ptr<jpcRes> _eta2JPC;

  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_EtaToPhiPhi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_F0ToPhiPhi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_F1ToPhiPhi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_F2ToPhiPhi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_Eta2ToPhiPhi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_PhiToKK;

};

