#ifndef _K1400Spin1Lh_H
#define _K1400Spin1Lh_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "Examples/Psi2SToKpKmPiGam/AbsPsi2SToKpKmPiGamLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class K1400Spin1Lh : public AbsPsi2SToKpKmPiGamLh{

public:

  // create/copy/destroy:

  ///Constructor 
  K1400Spin1Lh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList>, boost::shared_ptr<const Psi2SToKpKmPiGamStates>);
  K1400Spin1Lh(boost::shared_ptr<AbsPsi2SToKpKmPiGamLh>);

  /** Destructor */
  virtual ~K1400Spin1Lh();

  virtual AbsPsi2SToKpKmPiGamLh* clone_() const {
    return new K1400Spin1Lh(_Psi2SToKpKmPiGamEvtListPtr, _Psi2SToKpKmPiGamStatesPtr);
  }


  // Getters:
  
  virtual void setMnUsrParams(MnUserParameters& upar, Psi2SToKpKmPiGamData::fitParamVal& startVal,  Psi2SToKpKmPiGamData::fitParamVal& errVal);
  virtual int setFitParamVal(Psi2SToKpKmPiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const;

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(Psi2SToKpKmPiGamData::fitParamVal& theParamVal) const;

protected:


private:
  virtual  complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, const Psi2SToKpKmPiGamData::fitParamVal& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData);

};

#endif
