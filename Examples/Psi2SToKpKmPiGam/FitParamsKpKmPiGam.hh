#ifndef _FitParamsKpKmPiGam_H
#define _FitParamsKpKmPiGam_H

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
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

#include "Minuit2/MnUserParameters.h"

// using namespace std;
using namespace ROOT::Minuit2;

struct paramKpKmPiGam
{
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToChiGam;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK890K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_0_K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_1_K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_2_K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa0Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa2Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKappa_K;
  pair<double, double> BwK890;
  pair<double, double> BwK1400_0;
  pair<double, double> BwK1400_1;
  pair<double, double> BwK1400_2;
  pair<double, double> Bwa980;
  pair<double, double> Bwa2;
  pair<double, double> BwKappa;

  double FlatMa980;
  double FlatgKK;
  double FlatgEtaPi;
  double phaseSpace;
};

struct paramEnumKpKmPiGam{
  enum {ChiGam=0, K890K, K_0_1400K, K_1_1400K, K_2_1400K, a980Pi, KappaK, nAmps,
	K890=nAmps, K_0_1400, K_1_1400, K_2_1400, Kappa, nMasses,
	a980=nMasses,nFlatte,
	phaseSpace=nFlatte, nPhaseSpace};

  static const std::string& name(unsigned int t)
  {
    static std::string fitName[paramEnumKpKmPiGam::nPhaseSpace]
      ={"ChiGam", "K890K", "K_0_1400K", "K_1_1400K", "K_2_1400K", "a980Pi", "KappaK",
	"K890", "K_0_1400", "K_1_1400", "K_2_1400", "Kappa",
	"a980",
	"phaseSpace"};
    if (t<0 || t>=paramEnumKpKmPiGam::nPhaseSpace) assert(0);
    return fitName[t];
  } 
};



class FitParamsKpKmPiGam{

public:

  FitParamsKpKmPiGam(); 

  virtual ~FitParamsKpKmPiGam();

  std::vector< boost::shared_ptr<const JPCLS> >  jpclsVec(unsigned int index);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ampMap(paramKpKmPiGam& params, unsigned int index);
  pair<double, double>& massPair(paramKpKmPiGam& params, unsigned int index);

  int setFitParamValDec(paramKpKmPiGam& theParamVal, const std::vector<double>& par, int counter, unsigned int index);
  int setFitParamValMass(paramKpKmPiGam& theParamVal, const std::vector<double>& par, int counter, unsigned int index);
  int setFitParamFlattea980Mass(paramKpKmPiGam& theParamVal, const std::vector<double>& par, int counter, std::string key);
  void setMnUsrParamsDec(MnUserParameters& upar, paramKpKmPiGam& startVal,  paramKpKmPiGam& errVal, unsigned int index);
  void setMnUsrParamsMass(MnUserParameters& upar, paramKpKmPiGam& startVal,  paramKpKmPiGam& errVal, unsigned int index);
  void setMnUsrParamsFlattea980Mass(MnUserParameters& upar, paramKpKmPiGam& startVal,  paramKpKmPiGam& errVal, std::string key);

protected:

private:
  std::map < unsigned int, std::vector< boost::shared_ptr<const JPCLS> > > _jpclsMap;
};


#endif
