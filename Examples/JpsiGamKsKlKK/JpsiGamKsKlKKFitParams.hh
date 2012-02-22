#ifndef _JpsiGamKsKlKKFitParams_H
#define _JpsiGamKsKlKKFitParams_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <complex>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"


#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/FitParamsBase.hh"

#include "Utils/PawianCollectionUtils.hh"

#include "Minuit2/MnUserParameters.h"




// using namespace std;
using namespace ROOT::Minuit2;


struct paramEnumJpsiGamKsKlKK{
  enum { PsiToEtacGamma=0, PsiToEta2225Gamma, PsiToF02020Gamma, PsiToF22010Gamma, PsiToF22300Gamma,  PsiToF22340Gamma, PsiToEta21870Gamma, PsiToF1Gamma, 
	 F02020ToPhiPhi, F22300ToPhiPhi, Eta21870ToPhiPhi, F1ToPhiPhi, nAmps,
	 etac=nAmps, eta2225, f02020, f22010, f22300 ,f22340, eta21870, f1, nMasses, f02020gKK=nMasses, f02020gPhiPhi, ngFactors,
        phaseSpace=ngFactors, nOthers };
  
  static const std::string& name(unsigned int t)
  {
    static std::string fitName[paramEnumJpsiGamKsKlKK::nOthers]
      ={"PsiToEtacGam", "PsiToEta2225Gamma", "PsiToF02020Gamma", "PsiToF22010Gamma", "PsiToF22300Gamma", "PsiToF2340Gamma", "PsiToEta21870Gamma", "PsiToF1Gamma",
	"F02020ToPhiPhi", "F22300ToPhiPhi", "Eta21870ToPhiPhi", "F1ToPhiPhi",         
	"etac", "eta2250", "f02020", "f22010", "f22300", "f22340", "eta21870", "f1", "f02020gKK","f02020gPhiPhi", "phaseSpace"  };
    if (t<0 || t>=paramEnumJpsiGamKsKlKK::nOthers) assert(0);
    return fitName[t];
  } 
};



class JpsiGamKsKlKKFitParams : public FitParamsBase {

public:

  JpsiGamKsKlKKFitParams(); 
  JpsiGamKsKlKKFitParams(fitParams& theStartparams, fitParams& theErrorparams);
  virtual ~JpsiGamKsKlKKFitParams();
  
  virtual const std::string ampName(int index) {return paramEnumJpsiGamKsKlKK::name(index);}
  virtual const std::string massName(int index) {return paramEnumJpsiGamKsKlKK::name(index);}
  virtual const std::string widthName(int index) {return paramEnumJpsiGamKsKlKK::name(index);}
  virtual const std::string gFactorName(int index) {return paramEnumJpsiGamKsKlKK::name(index);}
  virtual const std::string otherName(int index) {return paramEnumJpsiGamKsKlKK::name(index);}

  virtual int ampIdxMin() {return paramEnumJpsiGamKsKlKK::PsiToEtacGamma;}
  virtual int ampIdxMax() {return paramEnumJpsiGamKsKlKK::nAmps-1;}
  virtual int massIdxMin() {return paramEnumJpsiGamKsKlKK::nAmps;}
  virtual int massIdxMax() {return paramEnumJpsiGamKsKlKK::nMasses-1;}
  virtual int gFactorIdxMin() {return paramEnumJpsiGamKsKlKK::nMasses;}
  virtual int gFactorIdxMax() {return paramEnumJpsiGamKsKlKK::ngFactors-1;}
  virtual int otherIdxMin() {return paramEnumJpsiGamKsKlKK::ngFactors;}
  virtual int otherIdxMax() {return paramEnumJpsiGamKsKlKK::nOthers-1;}

protected:

private:

};


#endif
