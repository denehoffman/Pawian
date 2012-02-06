#ifndef _FitParams2K2PiGam_H
#define _FitParams2K2PiGam_H

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
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

#include "Minuit2/MnUserParameters.h"

// using namespace std;
using namespace ROOT::Minuit2;

struct param2K2PiGam
{
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToChiGam;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2K892;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2Kappa;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2K_2_1430;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2K_0_1430;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_1430_K_2_1430;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_1_1410_K_1_1410;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_1_1400K;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K1400ToK892Pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_1_1270_K;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_1270ToK892Pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_1270ToK_0_1430Pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_1_2400K;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_2400Tof980K;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_pif1710_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_kf1710_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980f980;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980f2200;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980f_2_2300;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof_2_2300sigma;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf980;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf1370;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf1710;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf2200;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_pif1370_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_kf1370_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_pif1500_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_kf1500_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f1710_pif1370_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f1710_kf1370_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_pif_2_1270_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_kf_2_1270_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_pif_2_1430_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_kf_2_1430_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_pif_2_1525_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_kf_2_1525_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_pif_2_1950_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f980_kf_2_1950_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f1500_pif_2_1525_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f1500_kf_2_1525_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f1710_pif_2_1430_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f1710_kf_2_1430_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f1710_pif_2_1950_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > f1710_kf_2_1950_pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_2400ToKf980;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_2400ToKf_0_1710;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_1430K_0_1950;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKappaK_0_1430;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKappaK_0_1950;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK892K1680;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK892K2300;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_0_1460ToK892Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_0_1460ToK_0_1430Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_0_1830ToK892Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_0_1830ToK_0_1430Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_1_1650K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_1650ToK892Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_1650ToK_0_1430Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi_2_1670Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_1670Tof_2_1270Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_1670ToK892K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi1800Pi0Tof980;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi1800Pi0Tof1370;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi1800Pi0ToKappa;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi1800Pi0ToK892K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi3000Pi0ToK892K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi3000Pi0ToK_0_1950K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi_2_2285Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_2285Tof1700Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_2285ToK892KPi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_2285ToK_0_1430K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_2285ToK_2_1430K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_2_1770K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_2_1770ToK_2_1430Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_1430KPi;

    pair<double, double> BwK892;
    pair<double, double> BwKappa;
    pair<double, double> BwK_2_1430;
    pair<double, double> BwK_0_1430;
    pair<double, double> BwK_1_1400;
    pair<double, double> BwK_1_1410;
    pair<double, double> BwK_1_1270;
    pair<double, double> Bwf1500;
    pair<double, double> Bwf1710;
    pair<double, double> Bwf2200;
    pair<double, double> BwSigma;
    pair<double, double> Bwf1360;
    pair<double, double> Bwf1370;
    pair<double, double> Bwf_2_1270;
    pair<double, double> Bwf_2_1430;
    pair<double, double> Bwf_2_1525;
    pair<double, double> Bwf_2_1950;
    pair<double, double> Bwf_2_2300;
    pair<double, double> BwK_0_2400;
    pair<double, double> BwK_1_2400;
    pair<double, double> BwK_2_1950;
    pair<double, double> BwK_0_1950;
    pair<double, double> BwK_1_1680;
    pair<double, double> BwK_1_2300;
    pair<double, double> BwK_0_1460;
    pair<double, double> BwK_0_1830;
    pair<double, double> BwK_1_1650;
    pair<double, double> BwPi_2_1670;
  pair<double, double> BwPi1800;
  pair<double, double> BwPi3000;
  pair<double, double> BwPi_2_2285;
  pair<double, double> BwK_2_1770;

    double Flatf980;
    double Flatf980gPiPi;
    double Flatf980gKK;

    double phaseSpace;
};

struct paramEnum2K2PiGam{
  enum {ChiGam=0, K892K892, KappaKappa, K_0_1430K_0_1430, K_2_1430K_2_1430, K_0_1430K_2_1430, K_1_1410K_1_1410,
	ChiToK_1_1400K, K_1_1400ToK892Pi, ChiToK_1_1270_K, K_1_1270ToK892Pi, K_1_1270ToK_0_1430Pi, 
	ChiToK_1_1650K, K_1_1650ToK892Pi, K_1_1650ToK_0_1430Pi, K_0_1430K_0_1950, KappaK_0_1950, ChiToKappaK_0_1430,
	ChiToK_1_2400K, K_1_2400Tof980K,
	f980_pif1710_k, f980_kf1710_pi, f980f980, f980f2200, ChiTof980f_2_2300, ChiTof_2_2300sigma, 
	f980_pif1370_k, f980_kf1370_pi, f980_pif1500_k, f980_kf1500_pi, f1710_pif1370_k, f1710_kf1370_pi,
	K_0_2400KToKf980, ChiToK_0_2400ToKf_0_1710, K892K_1_1680, K892K_1_2300, 
	sigmaf980, ChiToSigmaf1370, sigmaf1710, sigmaf2200, K_0_1460ToK892Pi, K_0_1460ToK_0_1430Pi, K_0_1830ToK892Pi, K_0_1830ToK_0_1430Pi, 
	f980_pif_2_1270_k, f980_kf_2_1270_pi, f980_pif_2_1430_k, f980_kf_2_1430_pi, f980_pif_2_1525_k, f980_kf_2_1525_pi, f980_pif_2_1950_k, f980_kf_2_1950_pi,
	f1500_pif_2_1525_k, f1500_kf_2_1525_pi, f1710_pif_2_1430_k, f1710_kf_2_1430_pi, f1710_pif_2_1950_k, f1710_kf_2_1950_pi, 
	ChiToPi_2_1670Pi, Pi_2_1670Tof_2_1270Pi, Pi_2_1670ToK892K, 
	ChiToPi1800Pi0Tof980, ChiToPi1800Pi0Tof1370, ChiToPi1800Pi0ToKappa, ChiToPi1800Pi0ToK892K, ChiToPi3000Pi0ToK892K, ChiToPi3000Pi0ToK_0_1950K, 
	ChiToPi_2_2285Pi,Pi_2_2285Tof1700Pi, Pi_2_2285ToK892KPi, Pi_2_2285ToK_0_1430K, Pi_2_2285ToK_2_1430K, 
	ChiToK_2_1770K, K_2_1770ToK_2_1430Pi, ChiToK_0_1430KPi, nAmps,         

	K892=nAmps,Kappa, K_0_1430, K_1_1400,  K_1_1410, K_2_1430, K_1_1270,K_1_1650, 
	f1500, f1710, f2200, sigma, f1360, f1370, K_0_2400, K_1_2400, K_0_1950, K_1_1680, K_1_2300, K_0_1460, K_0_1830, 
	f_2_1270, f_2_1430, f_2_1525, f_2_1950, f_2_2300, 
	Pi_2_1670, Pi1800, Pi3000, Pi_2_2285, K_2_1770, nMasses,

	f980=nMasses,f980gpipi, f980gKK, nFlatte,

	phaseSpace=nFlatte, nPhaseSpace};

  static const std::string& name(unsigned int t)
  {
    static std::string fitName[paramEnum2K2PiGam::nPhaseSpace]
      ={"ChiGam", "K892K892", "KappaKappa", "K_0_1430K_0_1430", "K_2_1430K_2_1430", "K_0_1430K_2_1430", "K_1_1410K_1_1410",
	"ChiToK_1_1400K", "K_1_1400ToK892Pi", "ChiToK_1_1270_K", "K_1_1270ToK892Pi", "K_1_1270ToK_0_1430Pi", 
	"ChiToK_1_1650K", "K_1_1650ToK892Pi", "K_1_1650ToK_0_1430Pi","K_0_1430K_0_1950", "KappaK_0_1950","ChiToKappaK_0_1430",
	"ChiToK_1_2400K", "K_1_2400Tof980K",
	"f980_pif1710_k", "f980_kf1710_pi", "f980f980", "f980f2200", "ChiTof980f_2_2300","ChiTof_2_2300sigma",
	"f980_pif1370_k", "f980_kf1370_pi", "f980_pif1500_k", "f980_kf1500_pi", "f1710_pif1370_k", "f1710_kf1370_pi",
	"K_0_2400KToKf980", "ChiToK_0_2400ToKf_0_1710", "K892K_1_1680", "K892K_1_2300", 
	"sigmaf980", "ChiToSigmaf1370", "sigmaf1710", "sigmaf2200", "K_0_1460ToK892Pi", "K_0_1460ToK_0_1430Pi","K_0_1830ToK892Pi", "K_0_1830ToK_0_1430Pi",
	"f980_pif_2_1270_k", "f980_kf_2_1270_pi","f980_pif_2_1430_k", "f980_kf_2_1430_pi","f980_pif_2_1525_k", "f980_kf_2_1525_pi", "f980_pif_2_1950_k", "f980_kf_2_1950_pi",  
	"f1500_pif_2_1525_k", "f1500_kf_2_1525_pi", "f1710_pif_2_1430_k", "f1710_kf_2_1430_pi", "f1710_pif_2_1950_k", "f1710_kf_2_1950_pi",
	"ChiToPi_2_1670Pi","Pi_2_1670Tof_2_1270Pi","Pi_2_1670ToK892K",
	"ChiToPi1800Pi0Tof980","ChiToPi1800Pi0Tof1370","ChiToPi1800Pi0ToKappa","ChiToPi1800Pi0ToK892K","ChiToPi3000Pi0ToK892K","ChiToPi3000Pi0ToK_0_1950K",
	"ChiToPi_2_2285Pi","Pi_2_2285Tof1700Pi","Pi_2_2285ToK892KPi","Pi_2_2285ToK_0_1430K","Pi_2_2285ToK_2_1430K", 
	"ChiToK_2_1770K","K_2_1770ToK_2_1430Pi","ChiToK_0_1430KPi",

	"K892", "Kappa", "K_0_1430", "K_1_1400", "K_1_1410", "K_2_1430", "K_1_1270", "K_1_1650", 
        "f1500", "f1710", "f2200", "sigma", "f1360", "f1370", "K_0_2400", "K_1_2400", "K_0_1950", "K_1_1680", "K_1_2300", "K_0_1460", "K_0_1830",
	"f_2_1270", "f_2_1430", "f_2_1525","f_2_1950","f_2_2300",
	"Pi_2_1670","Pi1800","Pi3000","Pi_2_2285","K_2_1770",

	"f980","f980gpipi", "f980gKK",
	
	"phaseSpace"};
    if (t<0 || t>=paramEnum2K2PiGam::nPhaseSpace) assert(0);
    return fitName[t];
  } 
};



class FitParams2K2PiGam{

public:

  FitParams2K2PiGam(); 

  virtual ~FitParams2K2PiGam();

  std::vector< boost::shared_ptr<const JPCLS> >  jpclsVec(unsigned int index);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ampMap(param2K2PiGam& params, unsigned int index);
  pair<double, double>& massPair(param2K2PiGam& params, unsigned int index);

  int setFitParamValDec(param2K2PiGam& theParamVal, const std::vector<double>& par, int counter, unsigned int index);
  int setFitParamValMass(param2K2PiGam& theParamVal, const std::vector<double>& par, int counter, unsigned int index);
  int setFitParamFlatteMass(param2K2PiGam& theParamVal, const std::vector<double>& par, int counter, std::string key);
  void setMnUsrParamsDec(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal, unsigned int index);
  void setMnUsrParamsMass(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal, unsigned int index);
  void setMnUsrParamsFlatteMass(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal, std::string key);

protected:

private:
  std::map < unsigned int, std::vector< boost::shared_ptr<const JPCLS> > > _jpclsMap;
};


#endif
