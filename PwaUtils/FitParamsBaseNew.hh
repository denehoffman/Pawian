#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <complex>

#include <boost/shared_ptr.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

#include "Minuit2/MnUserParameters.h"

using namespace ROOT::Minuit2;

typedef std::map<std::string, std::map< boost::shared_ptr<const JPCLS>, double, 
					pawian::Collection::SharedPtrLess > > mapStrJPCLS;
typedef std::map<std::string, std::map< boost::shared_ptr<const JPClamlam>, double, 
					pawian::Collection::SharedPtrLess > > mapStrJPCLamLam;
typedef std::map<std::string, double> mapStrDouble;

struct fitParamsNew {
  mapStrJPCLS Mags;
  mapStrJPCLS Phis;
  mapStrJPCLamLam MagLamLams;
  mapStrJPCLamLam PhiLamLams;
  mapStrDouble Masses;
  mapStrDouble Widths;
  mapStrDouble gFactors;
  mapStrDouble otherParams;
};

class FitParamsBaseNew {

public:
  FitParamsBaseNew(); 

  virtual ~FitParamsBaseNew();

  virtual void setMnUsrParams(MnUserParameters& upar, fitParamsNew& theValParams, fitParamsNew& theErrParams);
  virtual void getFitParamVal(const std::vector<double>& par, fitParamsNew& theParams);
  virtual void printParams(fitParamsNew& theParams);
  virtual void dumpParams(std::ostream& os, fitParamsNew& theVals,  fitParamsNew& theErrs);

protected:
  const std::string _magSuffix;
  const std::string _phiSuffix;
  const std::string _massSuffix;
  const std::string _widthSuffix;
  const std::string _gFactorSuffix;
  const std::string _otherParamSuffix;

private:
  virtual void setMnUsrParamsJPCLamLam(MnUserParameters& upar, 
				       mapStrJPCLamLam& startLamLamMagMap, 
				       mapStrJPCLamLam& errLamLamMagMap, 
				       const std::string& suffix);
  virtual void setMnUsrParamsJPCLS(MnUserParameters& upar, mapStrJPCLS& startMagMap, mapStrJPCLS& errMagMap, const std::string& suffix);
  virtual void setMnUsrParamsDouble(MnUserParameters& upar, mapStrDouble& startDoubleMap, mapStrDouble& errDoubleMap, const std::string& suffix);
  virtual void printJPCLamLamParams(fitParamsNew& theParams, mapStrJPCLamLam& lamLamMagMap, const std::string& suffix);
  virtual void printJPCLSParams(fitParamsNew& theParams, mapStrJPCLS& JPCLSMagMap, const std::string& suffix);
  virtual void printDoubleParams(fitParamsNew& theParams, mapStrDouble& doubleMap,  const std::string& suffix);
  virtual void dumpLamLamParams(std::ostream& os, mapStrJPCLamLam& valLamLamMap, mapStrJPCLamLam& errLamLamMap, const std::string& suffix);
  virtual void dumpJPCLSParams(std::ostream& os, mapStrJPCLS& valLSMap, mapStrJPCLS& errLSMap, const std::string& suffix);
  virtual void dumpDoubleParams(std::ostream& os, mapStrDouble& valDoubleMap, mapStrDouble& errDoubleMap, const std::string& suffix);
  virtual void getFitParamValJPCLamLam(const std::vector<double>& par, mapStrJPCLamLam& startLamLamMagMap, int& counter);
  virtual void getFitParamValJPCLS(const std::vector<double>& par, mapStrJPCLS& jpclsMap, int& counter);
  virtual void getFitParamValDouble(const std::vector<double>& par, mapStrDouble& doubleMap, int& counter);
};
