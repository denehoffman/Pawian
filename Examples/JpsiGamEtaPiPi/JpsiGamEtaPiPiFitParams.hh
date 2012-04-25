#ifndef _JpsiGamEtaPiPiFitParams_H
#define _JpsiGamEtaPiPiFitParams_H

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
#include "PwaUtils/FitParamsBase.hh"

// using namespace std;
using namespace ROOT::Minuit2;


struct paramEnumJpsiGamEtaPiPi{
  enum { PsiToEtaGamma=0, PsiToF1Gamma, nAmps,
	  	 a980=nAmps, nMasses, ngFactors=nMasses, phaseSpace=ngFactors, nOthers };
  
  static const std::string& name(unsigned int t)
  {
    static std::string fitName[paramEnumJpsiGamEtaPiPi::nOthers]
      ={"PsiToEtaGamma", "PsiToF1Gamma",
    		  "a980", "phaseSpace"};
    if (t<0 || t>=paramEnumJpsiGamEtaPiPi::nOthers) assert(0);
    return fitName[t];
  } 
};



class JpsiGamEtaPiPiFitParams : public FitParamsBase {

public:

  JpsiGamEtaPiPiFitParams();
  JpsiGamEtaPiPiFitParams(fitParams& theStartparams, fitParams& theErrorparams);
  virtual ~JpsiGamEtaPiPiFitParams();
  
  virtual const std::string ampName(int index) {return paramEnumJpsiGamEtaPiPi::name(index);}
  virtual const std::string massName(int index) {return paramEnumJpsiGamEtaPiPi::name(index);}
  virtual const std::string widthName(int index) {return paramEnumJpsiGamEtaPiPi::name(index);}
  virtual const std::string gFactorName(int index) {return paramEnumJpsiGamEtaPiPi::name(index);}
  virtual const std::string otherName(int index) {return paramEnumJpsiGamEtaPiPi::name(index);}

  virtual int ampIdxMin() {return paramEnumJpsiGamEtaPiPi::PsiToEtaGamma;}
  virtual int ampIdxMax() {return paramEnumJpsiGamEtaPiPi::nAmps-1;}
  virtual int massIdxMin() {return paramEnumJpsiGamEtaPiPi::nAmps;}
  virtual int massIdxMax() {return paramEnumJpsiGamEtaPiPi::nMasses-1;}
  virtual int gFactorIdxMin() {return paramEnumJpsiGamEtaPiPi::nMasses;}
  virtual int gFactorIdxMax() {return paramEnumJpsiGamEtaPiPi::ngFactors-1;}
  virtual int otherIdxMin() {return paramEnumJpsiGamEtaPiPi::ngFactors;}
  virtual int otherIdxMax() {return paramEnumJpsiGamEtaPiPi::nOthers-1;}

protected:

private:

};


#endif
