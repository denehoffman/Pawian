#ifndef _PsiToChic1GamProdLh_H
#define _PsiToChic1GamProdLh_H

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

#include "PwaUtils/AbsLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEvtList.hh"
#include "PwaUtils/DataUtils.hh"


class PsiToChic1GamProdLh : public AbsLh{

public:

  // create/copy/destroy:

  ///Constructor 
  PsiToChic1GamProdLh(boost::shared_ptr<const EvtDataBaseList>);
  PsiToChic1GamProdLh(boost::shared_ptr<AbsLh>);

  /** Destructor */
  virtual ~PsiToChic1GamProdLh();

  virtual AbsLh* clone_() const {
    return new PsiToChic1GamProdLh(_evtListPtr);
  }


  // Getters:
  virtual double calcEvtIntensity(EvtData* theData, fitParams& theParamVal);
  
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual void print(std::ostream& os) const;

protected:

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, fitParams& theParamVal, EvtData* theData);  
  complex<double> calcCoherentAmp(Spin Minit, Spin MChi, Spin MGamma, fitParams& theParamVal, EvtData* theData);

private:

};

#endif
