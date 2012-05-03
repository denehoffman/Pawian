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
#include "Examples/JpsiGamEtaPiPi/PsiToXGamStates.hh"


class JpsiGamEtaPiPiStates : public JpsiToXGamStates {

public:

  ///Constructor 
  JpsiGamEtaPiPiStates();


  /** Destructor */
  virtual ~JpsiGamEtaPiPiStates();


  std::vector< boost::shared_ptr<const JPCLS> > EtaToa0PiStates() const {return _JPCLS_EtaToa0Pi;}
  std::vector< boost::shared_ptr<const JPCLS> > EtaTof0EtaStates() const {return _JPCLS_EtaTof0Eta;}
  std::vector< boost::shared_ptr<const JPCLS> > EtaToa2PiStates() const {return _JPCLS_EtaToa2Pi;}
  std::vector< boost::shared_ptr<const JPCLS> > Eta2Toa0PiStates() const {return _JPCLS_Eta2Toa0Pi;}
  std::vector< boost::shared_ptr<const JPCLS> > Eta2Toa2PiStates() const {return _JPCLS_Eta2Toa2Pi;}
  std::vector< boost::shared_ptr<const JPCLS> > Eta2Tof0EtaStates() const {return _JPCLS_Eta2Tof0Eta;}
  std::vector< boost::shared_ptr<const JPCLS> > F1Toa0PiStates() const {return _JPCLS_F1Toa0Pi;}
  std::vector< boost::shared_ptr<const JPCLS> > F1Tof0EtaStates() const {return _JPCLS_F1Tof0Eta;} 
  std::vector< boost::shared_ptr<const JPCLS> > F1Toa2PiStates() const {return _JPCLS_F1Toa2Pi;}
 
  virtual void print(std::ostream& os) const; 

protected:


private:

  boost::shared_ptr<jpcRes> _etaJPC;
  boost::shared_ptr<jpcRes> _eta2JPC;
  boost::shared_ptr<jpcRes> _f0JPC;
  boost::shared_ptr<jpcRes> _f1JPC;
  boost::shared_ptr<jpcRes> _f2JPC;
  boost::shared_ptr<jpcRes> _a0JPC;
  boost::shared_ptr<jpcRes> _a2JPC;
  boost::shared_ptr<jpcRes> _piJPC;

  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_EtaToa0Pi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_EtaTof0Eta;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_F1Toa0Pi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_EtaToa2Pi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_Eta2Toa0Pi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_Eta2Toa2Pi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_Eta2Tof0Eta;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_F1Tof0Eta;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_F1Toa2Pi;  
};

#endif
