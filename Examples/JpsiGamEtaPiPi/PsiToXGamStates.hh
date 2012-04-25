#ifndef _JpsiToXGamStates_H
#define _JpsiToXGamStates_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>


#include <boost/shared_ptr.hpp>

// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"


class JpsiToXGamStates {

public:

  ///Constructor 
  JpsiToXGamStates();


  /** Destructor */
  virtual ~JpsiToXGamStates();

  virtual std::vector< boost::shared_ptr<const JPCLS> > PsiToEtaGammaStates() const {return _JPCLS_PsiToEtaGamma;}
  virtual std::vector< boost::shared_ptr<const JPCLS> > PsiToF0GammaStates() const {return _JPCLS_PsiToF0Gamma;}
  virtual std::vector< boost::shared_ptr<const JPCLS> > PsiToF1GammaStates() const {return _JPCLS_PsiToF1Gamma;}
  virtual std::vector< boost::shared_ptr<const JPCLS> > PsiToF2GammaStates() const {return _JPCLS_PsiToF2Gamma;}
  virtual std::vector< boost::shared_ptr<const JPCLS> > PsiToEta2GammaStates() const {return _JPCLS_PsiToEta2Gamma;}

  
  virtual void print(std::ostream& os) const; 

protected:
  virtual void fillJPCLS(boost::shared_ptr<jpcRes>, boost::shared_ptr<jpcRes>, boost::shared_ptr<jpcRes>, std::vector< boost::shared_ptr<const JPCLS> >& );
  virtual void printDecayJPCLS(std::ostream& os,std::vector< boost::shared_ptr<const JPCLS > > theJPCLS) const;

private:

  boost::shared_ptr<jpcRes> _psiJPC;
  boost::shared_ptr<jpcRes> _gammaJPC;
  boost::shared_ptr<jpcRes> _etaJPC;
  boost::shared_ptr<jpcRes> _f0JPC;
  boost::shared_ptr<jpcRes> _f1JPC;
  boost::shared_ptr<jpcRes> _f2JPC;
  boost::shared_ptr<jpcRes> _eta2JPC;
  
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_PsiToEtaGamma;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_PsiToF0Gamma;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_PsiToF1Gamma;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_PsiToF2Gamma;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_PsiToEta2Gamma;
  
};

#endif
