#ifndef _Psi2SToKpKmPiGamHyp1Lh_H
#define _Psi2SToKpKmPiGamHyp1Lh_H

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

#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamBaseLh.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class Psi2SToKpKmPiGamHyp1Lh : public Psi2SToKpKmPiGamBaseLh{

public:

  // create/copy/destroy:

  ///Constructor 
  Psi2SToKpKmPiGamHyp1Lh(boost::shared_ptr<const EvtDataBaseList>, const std::map<const std::string, bool>& hypMap);
  Psi2SToKpKmPiGamHyp1Lh(boost::shared_ptr<AbsLh>, const std::map<const std::string, bool>& hypMap);

  /** Destructor */
  virtual ~Psi2SToKpKmPiGamHyp1Lh();

  virtual AbsLh* clone_() const {
    return new Psi2SToKpKmPiGamHyp1Lh(_evtListPtr, _hypMap);
  }


  // Getters:
  virtual void print(std::ostream& os) const;

protected:

  complex<double> calcDecAmp(complex<double>& inAmp,Spin lamChi, fitParams& theParamVal, EvtData* theData);

private:
};

#endif
