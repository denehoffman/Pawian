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
  enum { PsiToEtaGamma=0, PsiToEta2Gamma, PsiToF1Gamma, 
	 EtaToPiPiEta, EtaToA980Pi, EtaToA2_1320Pi, EtaToSigmaEta, EtaTof0_980Eta, 
	 Eta2ToPiPiEta, Eta2ToA980Pi, Eta2ToA2_1320Pi, Eta2ToSigmaEta, Eta2Tof0_980Eta,
	 F1ToPiPiEta, F1ToA980Pi, F1ToA2_1320Pi, F1ToSigmaEta, F1Tof0_980Eta, nAmps,
	 a0_980=nAmps, a2_1320, sigma, f0_980, nMasses, 
	 a0_980gPiEta=nMasses, a0_980gKK, f0_980gPiPi, f0_980gKK, ngFactors, 
	 phaseSpace=ngFactors, nOthers };
  
  static const std::string& name(unsigned int t)
  {
    static std::string fitName[paramEnumJpsiGamEtaPiPi::nOthers]
      ={"PsiToEtaGamma", "PsiToEta2Gamma", "PsiToF1Gamma", 
	"EtaToPiPiEta", "EtaToA980Pi", "EtaToA2_1320Pi", "EtaToSigmaEta", "EtaTof0_980Eta", 
	"Eta2ToPiPiEta", "Eta2ToA980Pi", "Eta2ToA2_1320Pi", "Eta2ToSigmaEta", "Eta2Tof0_980Eta",
	"F1ToPiPiEta","F1ToA980Pi", "F1ToA2_1320Pi", "F1ToSigmaEta", "F1Tof0_980Eta",
	"a0_980","a2_1320","sigma","f0_980", 
	"a0_980gPiEta", "a0_980gKK","f0_980gPiPi", "f0_980gKK", 
	"phaseSpace"};
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
  virtual const std::string ampLamLamName(int index) {return paramEnumJpsiGamEtaPiPi::name(index);}
  virtual const std::string massName(int index) {return paramEnumJpsiGamEtaPiPi::name(index);}
  virtual const std::string widthName(int index) {return paramEnumJpsiGamEtaPiPi::name(index);}
  virtual const std::string gFactorName(int index) {return paramEnumJpsiGamEtaPiPi::name(index);}
  virtual const std::string otherName(int index) {return paramEnumJpsiGamEtaPiPi::name(index);}

  virtual int ampLamLamIdxMin() {return paramEnumJpsiGamEtaPiPi::PsiToEtaGamma;}
  virtual int ampLamLamIdxMax() {return paramEnumJpsiGamEtaPiPi::PsiToF1Gamma;}
  virtual int ampIdxMin() {return paramEnumJpsiGamEtaPiPi::EtaToPiPiEta;}
  virtual int ampIdxMax() {return paramEnumJpsiGamEtaPiPi::nAmps-1;}
  virtual int massIdxMin() {return paramEnumJpsiGamEtaPiPi::nAmps;}
  virtual int massIdxMax() {return paramEnumJpsiGamEtaPiPi::nMasses-1;}
  virtual int gFactorIdxMin() {return paramEnumJpsiGamEtaPiPi::nMasses;}
  virtual int gFactorIdxMax() {return paramEnumJpsiGamEtaPiPi::ngFactors-1;}
  virtual int otherIdxMin() {return paramEnumJpsiGamEtaPiPi::ngFactors;}
  virtual int otherIdxMax() {return paramEnumJpsiGamEtaPiPi::nOthers-1;}

protected:
  virtual void initialize();
private:

};


#endif
