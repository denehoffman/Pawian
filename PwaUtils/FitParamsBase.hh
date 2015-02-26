//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <complex>

#include <memory>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

#include "FitParams/AbsPawianParameters.hh"

//using namespace ROOT::Minuit2;

typedef std::map<std::string, std::map< std::shared_ptr<const LScomb>, double,
					pawian::Collection::SharedPtrLess > > mapStrLS;
typedef std::map<std::string, std::map< std::shared_ptr<const jpcRes>, double,
					pawian::Collection::SharedPtrLess > > mapStrJPC;
typedef std::map<std::string, std::map< std::shared_ptr<const JPCLS>, double,
					pawian::Collection::SharedPtrLess > > mapStrJPCLS;
typedef std::map<std::string, std::map< std::shared_ptr<const JPClamlam>, double,
					pawian::Collection::SharedPtrLess > > mapStrJPCLamLam;
typedef std::map<std::string, double> mapStrDouble;

struct fitParams {
  mapStrJPC MagsJPC;
  mapStrJPC PhisJPC;
  mapStrLS MagsLS;
  mapStrLS PhisLS;
  mapStrJPCLS Mags;
  mapStrJPCLS Phis;
  mapStrJPCLamLam MagLamLams;
  mapStrJPCLamLam PhiLamLams;
  mapStrDouble Masses;
  mapStrDouble Widths;
  mapStrDouble gFactors;
  mapStrDouble otherParams;
};

class FitParamsBase {

public:
  FitParamsBase();

  virtual ~FitParamsBase();

  virtual void setAbsPawianParams(std::shared_ptr<AbsPawianParameters> upar, fitParams& theValParams, fitParams& theErrParams);
  virtual void getFitParamVal(const std::vector<double>& par, fitParams& theParams);
  virtual void printParams(fitParams& theParams);
  virtual void dumpParams(std::ostream& os, fitParams& theVals,  fitParams& theErrs);

protected:
  //  const std::string _isoSuffix;
  const std::string _magSuffix;
  const std::string _phiSuffix;
  const std::string _massSuffix;
  const std::string _widthSuffix;
  const std::string _gFactorSuffix;
  const std::string _otherParamSuffix;

private:
  virtual void setAbsPawianParamsJPC(std::shared_ptr<AbsPawianParameters> upar,
				     mapStrJPC& startJPCMap,
				     mapStrJPC& errJPCMap,
				     const std::string& suffix);
  virtual void setAbsPawianParamsLS(std::shared_ptr<AbsPawianParameters> upar, mapStrLS& startMagMap, mapStrLS& errMagMap, const std::string& suffix);
  virtual void setAbsPawianParamsJPCLamLam(std::shared_ptr<AbsPawianParameters> upar,
				       mapStrJPCLamLam& startLamLamMagMap,
				       mapStrJPCLamLam& errLamLamMagMap,
				       const std::string& suffix);
  virtual void setAbsPawianParamsJPCLS(std::shared_ptr<AbsPawianParameters> upar, mapStrJPCLS& startMagMap, mapStrJPCLS& errMagMap, const std::string& suffix);
  virtual void setAbsPawianParamsDouble(std::shared_ptr<AbsPawianParameters> upar, mapStrDouble& startDoubleMap, mapStrDouble& errDoubleMap, const std::string& suffix);
  virtual void printJPCParams(fitParams& theParams, mapStrJPC& jpcMap, const std::string& suffix);
  virtual void printLSParams(fitParams& theParams, mapStrLS& LSMagMap, const std::string& suffix);
  virtual void printJPCLamLamParams(fitParams& theParams, mapStrJPCLamLam& lamLamMagMap, const std::string& suffix);
  virtual void printJPCLSParams(fitParams& theParams, mapStrJPCLS& JPCLSMagMap, const std::string& suffix);
  virtual void printDoubleParams(fitParams& theParams, mapStrDouble& doubleMap,  const std::string& suffix);
  virtual void dumpJPCParams(std::ostream& os, mapStrJPC& valJPCMap, mapStrJPC& errJPCMap, const std::string& suffix);
  virtual void dumpLSParams(std::ostream& os, mapStrLS& valLSMap, mapStrLS& errLSMap, const std::string& suffix);
  virtual void dumpLamLamParams(std::ostream& os, mapStrJPCLamLam& valLamLamMap, mapStrJPCLamLam& errLamLamMap, const std::string& suffix);
  virtual void dumpJPCLSParams(std::ostream& os, mapStrJPCLS& valLSMap, mapStrJPCLS& errLSMap, const std::string& suffix);
  virtual void dumpDoubleParams(std::ostream& os, mapStrDouble& valDoubleMap, mapStrDouble& errDoubleMap, const std::string& suffix);
  virtual void getFitParamValJPC(const std::vector<double>& par, mapStrJPC& startMap, unsigned int& counter);
  virtual void getFitParamValLS(const std::vector<double>& par, mapStrLS& lsMap, unsigned int& counter);
  virtual void getFitParamValJPCLamLam(const std::vector<double>& par, mapStrJPCLamLam& startLamLamMagMap, unsigned int& counter);
  virtual void getFitParamValJPCLS(const std::vector<double>& par, mapStrJPCLS& jpclsMap, unsigned int& counter);
  virtual void getFitParamValDouble(const std::vector<double>& par, mapStrDouble& doubleMap, unsigned int& counter);
};
