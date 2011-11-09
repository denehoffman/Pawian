#ifndef _FitParamsChic1ToKpKmPiGam_H
#define _FitParamsChic1ToKpKmPiGam_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <complex>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
// #include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "Utils/PawianCollectionUtils.hh"

// #include "Minuit2/MnUserParameters.h"

// using namespace std;
// using namespace ROOT::Minuit2;

struct paramEnumChic1ToKpKmPiGam{
  enum {ChiGam=0, K890K, K_0_1400K, K_1_1400K, K_2_1400K, K_1_1680K, a980Pi, KappaK, KKPi, nAmps,
	K890=nAmps, K_0_1400, K_1_1400, K_2_1400,  K_1_1680, Kappa, a980, nMasses,
	a980gKK=nMasses, a980gEtaPi, ngFactors,
	phaseSpace=ngFactors, nOthers};

  static const std::string& name(unsigned int t)
  {
    static std::string fitName[paramEnumChic1ToKpKmPiGam::nOthers]
      ={"ChiGam", "K890K", "K_0_1400K", "K_1_1400K", "K_2_1400K",  "K_1_1680K", "a980Pi", "KappaK", "KKPi",
	"K890", "K_0_1400", "K_1_1400", "K_2_1400", "K_1_1680", "Kappa",
	"a980",
	"a980gKK", "a980gEtaPi",
	"phaseSpace"};
    if (t<0 || t>=paramEnumChic1ToKpKmPiGam::nOthers) assert(0);
    return fitName[t];
  }

//   static int theenum(std::string theString){
//     if (theString=="K_1_1400") return paramEnumChic1ToKpKmPiGam::K_1_1400;
//     return -1;
//   }
};



class FitParamsChic1ToKpKmPiGam : public FitParamsBase{

public:

  FitParamsChic1ToKpKmPiGam(); 
  FitParamsChic1ToKpKmPiGam(fitParams& theStartparams, fitParams& theErrorparams); 

  virtual ~FitParamsChic1ToKpKmPiGam();

  virtual const std::string ampName(int index) {return paramEnumChic1ToKpKmPiGam::name(index);}
  virtual const std::string massName(int index) {return paramEnumChic1ToKpKmPiGam::name(index);}
  virtual const std::string widthName(int index) {return paramEnumChic1ToKpKmPiGam::name(index);}
  virtual const std::string gFactorName(int index) {return paramEnumChic1ToKpKmPiGam::name(index);}
  virtual const std::string otherName(int index) {return paramEnumChic1ToKpKmPiGam::name(index);}

protected:

private:

};


#endif
