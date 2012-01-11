#ifndef _Psi2STo2K2PiGamStates_H
#define _Psi2STo2K2PiGamStates_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>


#include <boost/shared_ptr.hpp>

// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"


class Psi2STo2K2PiGamStates {

public:

  ///Constructor 
  Psi2STo2K2PiGamStates();


  /** Destructor */
  virtual ~Psi2STo2K2PiGamStates();

  std::vector< boost::shared_ptr<const JPCLS> > PsiToChiGamStates() const {return _PsiJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiToK1400ToK892piStates() const {return _ChiToK1400ToK892piJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiTo2K892States() const {return _ChiTo2K892JPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiTo2K_2_1430States() const {return _ChiTo2K_2_1430JPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiTo2K_0_States() const {return _ChiTo2K_0_JPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiToK0K2_States() const {return _ChiToK0K2_JPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiTof0f0States() const {return _ChiTof0f0JPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiTof0f2States() const {return _ChiTof0f2JPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > K1400ToKst1PiStates() const {return _K1400ToKst1PiJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiToK0K0States() const {return _ChiToK0K0JPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > K1ToK0PiStates() const {return _K1ToK0PiJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiToPi_2PiStates() const {return _ChiToPi_2PiJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > Pi_2Tof_2PiStates() const {return _Pi_2Tof_2PiJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > Pi_2ToKst1KStates() const {return _Pi_2ToKst1KJPCLS;}
  std::vector< boost::shared_ptr<const JPCLS> > ChiToPi0Pi0States() const {return _ChiToPi0Pi0JPCLS;}
  void print(std::ostream& os) const; 

protected:


private:
  boost::shared_ptr<jpcRes> _psiJPC;
  boost::shared_ptr<jpcRes> _chic0JPC;
  boost::shared_ptr<jpcRes> _Kst0JPC;
  boost::shared_ptr<jpcRes> _Kst1JPC;
  boost::shared_ptr<jpcRes> _Kst2JPC;
  boost::shared_ptr<jpcRes> _K1400JPC;
  boost::shared_ptr<jpcRes> _K0plusJPC;
  boost::shared_ptr<jpcRes> _f0JPC;
  boost::shared_ptr<jpcRes> _f2JPC;
  boost::shared_ptr<jpcRes> _pi0JPC;
  boost::shared_ptr<jpcRes> _kJPC;
  boost::shared_ptr<jpcRes> _gamJPC;
  boost::shared_ptr<jpcRes> _pi2JPC;

  std::vector< boost::shared_ptr<const JPCLS> > _PsiJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiToK1400ToK892piJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiTo2K892JPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiTo2K_2_1430JPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiTo2K_0_JPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiToK0K2_JPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiTof0f0JPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiTof0f2JPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiToK0K0JPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _K1ToK0PiJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> >  _K1400ToKst1PiJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiToPi_2PiJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _Pi_2Tof_2PiJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _Pi_2ToKst1KJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _ChiToPi0Pi0JPCLS;

  std::vector< boost::shared_ptr<const JPCLS> > _Kst1JPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _Kst2JPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > _f2JPCLS;

  void fillJPCLS(boost::shared_ptr<jpcRes>, boost::shared_ptr<jpcRes>, boost::shared_ptr<jpcRes>, std::vector< boost::shared_ptr<const JPCLS> >& ); 
};

#endif
