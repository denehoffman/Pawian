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

#include "Examples/Psi2SToKpKmPiGam/AbsPsi2SToKpKmPiGamLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class PsiToChic1GamProdLh : public AbsPsi2SToKpKmPiGamLh{

public:

  // create/copy/destroy:

  ///Constructor 
  PsiToChic1GamProdLh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList>, boost::shared_ptr<const Psi2SToKpKmPiGamStates>);
  PsiToChic1GamProdLh(boost::shared_ptr<AbsPsi2SToKpKmPiGamLh>);

  /** Destructor */
  virtual ~PsiToChic1GamProdLh();

  virtual AbsPsi2SToKpKmPiGamLh* clone_() const {
    return new PsiToChic1GamProdLh(_Psi2SToKpKmPiGamEvtListPtr, _Psi2SToKpKmPiGamStatesPtr);
  }


  // Getters:
  virtual double calcEvtIntensity(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, const Psi2SToKpKmPiGamData::fitParamVal& theParamVal);
  
  virtual void setMnUsrParams(MnUserParameters& upar, Psi2SToKpKmPiGamData::fitParamVal& startVal,  Psi2SToKpKmPiGamData::fitParamVal& errVal);
  virtual int setFitParamVal(Psi2SToKpKmPiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const;

  virtual void print(std::ostream& os) const;
  //  virtual void printCurrentFitResult(Psi2SToKpKmPiGamData::fitParamVal& theParamVal) const;

protected:


private:

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, const Psi2SToKpKmPiGamData::fitParamVal& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData);  
  complex<double> calcCoherentAmp(Spin Minit, Spin MChi, Spin MGamma, const Psi2SToKpKmPiGamData::fitParamVal& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData);

};

#endif
