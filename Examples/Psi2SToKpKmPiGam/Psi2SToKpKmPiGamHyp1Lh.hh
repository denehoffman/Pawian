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
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class Psi2SToKpKmPiGamHyp1Lh : public Psi2SToKpKmPiGamBaseLh{

public:

  // create/copy/destroy:

  ///Constructor 
  Psi2SToKpKmPiGamHyp1Lh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList>, const std::map<const std::string, bool>& hypMap);
  Psi2SToKpKmPiGamHyp1Lh(boost::shared_ptr<AbsPsi2SToKpKmPiGamLh>, const std::map<const std::string, bool>& hypMap);

  /** Destructor */
  virtual ~Psi2SToKpKmPiGamHyp1Lh();

  virtual AbsPsi2SToKpKmPiGamLh* clone_() const {
    return new Psi2SToKpKmPiGamHyp1Lh(_Psi2SToKpKmPiGamEvtListPtr, _hypMap);
  }


  // Getters:
  
  virtual void setMnUsrParams(MnUserParameters& upar, paramKpKmPiGam& startVal,  paramKpKmPiGam& errVal);
  virtual int setFitParamVal(paramKpKmPiGam& theParamVal, const std::vector<double>& par);

  virtual unsigned int nFitParams();
  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(paramKpKmPiGam& theParamVal);
  virtual void dumpCurrentResult(std::ostream& os, paramKpKmPiGam& theParamVal, std::string& suffix);
protected:

  complex<double> calcDecAmp(complex<double>& inAmp,Spin lamChi, const paramKpKmPiGam& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData);

private:
  unsigned int _nFitParams;
  std::vector<unsigned int> _ampVec;
  std::vector<unsigned int> _massVec;


};

#endif
