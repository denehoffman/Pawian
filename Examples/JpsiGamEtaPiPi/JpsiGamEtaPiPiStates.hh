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

  
  virtual void print(std::ostream& os) const; 

protected:


private:

  boost::shared_ptr<jpcRes> _etaJPC;
  boost::shared_ptr<jpcRes> _f0JPC;
  boost::shared_ptr<jpcRes> _f1JPC;
  boost::shared_ptr<jpcRes> _f2JPC;
  boost::shared_ptr<jpcRes> _a0JPC;
  boost::shared_ptr<jpcRes> _piJPC;

  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_EtaToa0Pi;

  
};

#endif
