#ifndef _FitParamIndex_H
#define _FitParamIndex_H

#include "PwaUtils/FitParamsBase.hh"

class FitParamIndex {
  
  public:
  
  FitParamIndex( fitParams& theParams  );
  virtual ~FitParamIndex(){};
  
  int Mag(int amp, boost::shared_ptr<const JPCLS> state){ return _Mags[amp][state]; }
  int Phi(int amp, boost::shared_ptr<const JPCLS> state){ return _Phis[amp][state]; }
  
  
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, int, pawian::Collection::SharedPtrLess > > _Mags;
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, int, pawian::Collection::SharedPtrLess > > _Phis;
  

};

#endif
