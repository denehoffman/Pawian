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

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKData.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

#include "Minuit2/MnUserParameters.h"

// using namespace std;
using namespace ROOT::Minuit2;

struct paramGamKsKlKK
{
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToEtacGamma;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToEta2225Gamma;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToF02020Gamma;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > F02020ToPhiPhi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToF22010Gamma;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToF22300Gamma;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToF22340Gamma;

  
  pair<double, double> BwEtac;
  pair<double, double> BwEta2225;
  pair<double, double> BwF02020;
  pair<double, double> BwF22010;
  pair<double, double> BwF22300;
  pair<double, double> BwF22340;

};

struct paramEnumGamKsKlKK{
  enum { etacGamma=0, eta2225Gamma, f02020Gamma, f22010Gamma, f22300Gamma,  f22340Gamma, 
	 f02020ToPhiPhi, nAmps,
	 etac=nAmps, eta2225, f02020, f22010, f22300 ,f22340, nMasses };
  
  static const std::string& name(unsigned int t)
  {
    static std::string fitName[paramEnumGamKsKlKK::nMasses]
      ={"etacGam", "eta2225Gamma", "f02020Gamma", "f22010Gamma", "f22300Gamma", "f2340Gamma",
	"f02020ToPhiPhi",           
	"etac", "eta2250", "f02020", "f22010", "f22340", "f22340" };
    if (t<0 || t>=paramEnumGamKsKlKK::nMasses) assert(0);
    return fitName[t];
  } 
};



class JpsiGamKsKlKKFitParams{

public:

  JpsiGamKsKlKKFitParams(); 

  virtual ~JpsiGamKsKlKKFitParams();

  std::vector< boost::shared_ptr<const JPCLS> >  jpclsVec(unsigned int index);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ampMap(paramGamKsKlKK& params, unsigned int index);
  pair<double, double>& massPair(paramGamKsKlKK& params, unsigned int index);

  int setFitParamValDec(paramGamKsKlKK& theParamVal, const std::vector<double>& par, int counter, unsigned int index);
  int setFitParamValMass(paramGamKsKlKK& theParamVal, const std::vector<double>& par, int counter, unsigned int index);
  void setMnUsrParamsDec(MnUserParameters& upar, paramGamKsKlKK& startVal,  paramGamKsKlKK& errVal, unsigned int index);
  void setMnUsrParamsMass(MnUserParameters& upar, paramGamKsKlKK& startVal,  paramGamKsKlKK& errVal, unsigned int index);
  
protected:

private:
  std::map < unsigned int, std::vector< boost::shared_ptr<const JPCLS> > > _jpclsMap;
};


#endif
