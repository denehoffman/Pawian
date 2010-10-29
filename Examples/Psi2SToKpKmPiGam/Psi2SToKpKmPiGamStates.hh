#ifndef _Psi2SToKpKmPiGamStates_H
#define _Psi2SToKpKmPiGamStates_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>


#include <boost/shared_ptr.hpp>

// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"


class Psi2SToKpKmPiGamStates {

public:

  ///Constructor 
  Psi2SToKpKmPiGamStates();


  /** Destructor */
  virtual ~Psi2SToKpKmPiGamStates();

  std::vector< boost::shared_ptr<const JPCLS> > PsiToChiGamStates() const {return _PsiJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiToKst0KStates() const {return _ChiToKst0KJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiToKst1KStates() const {return _ChiToKst1KJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiToKst2KStates() const {return _ChiToKst2KJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiToa0PiStates() const {return _ChiToa0PiJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiToa2PiStates() const {return _ChiToa2PiJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > Kst1ToKPiStates() const {return _Kst1JPCLS;} 
  std::vector< boost::shared_ptr<const JPCLS> > Kst2ToKPiStates() const {return _Kst2JPCLS;}  
  std::vector< boost::shared_ptr<const JPCLS> > a980ToKKStates() const {return _a0JPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > a2ToKKStates() const {return _a2JPCLS;}  

  void print(std::ostream& os) const; 

protected:


private:
  boost::shared_ptr<jpcRes> _psiJPC;
  boost::shared_ptr<jpcRes> _chic1JPC;
  boost::shared_ptr<jpcRes> _Kst0JPC;
  boost::shared_ptr<jpcRes> _Kst1JPC;
  boost::shared_ptr<jpcRes> _Kst2JPC;
  boost::shared_ptr<jpcRes> _a0JPC;
  boost::shared_ptr<jpcRes> _a2JPC;
  boost::shared_ptr<jpcRes> _pi0JPC;
  boost::shared_ptr<jpcRes> _kJPC;
  boost::shared_ptr<jpcRes> _gamJPC;

  std::vector< boost::shared_ptr<const JPCLS> > _PsiJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiToKst0KJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiToKst1KJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiToKst2KJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiToa0PiJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiToa2PiJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _Kst1JPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _Kst2JPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _a0JPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _a2JPCLS;

  void fillJPCLS(boost::shared_ptr<jpcRes>, boost::shared_ptr<jpcRes>, boost::shared_ptr<jpcRes>, std::vector< boost::shared_ptr<const JPCLS> >& ); 
};

#endif
