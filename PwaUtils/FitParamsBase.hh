#ifndef _FitParamsBase_H
#define _FitParamsBase_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <complex>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

#include "Minuit2/MnUserParameters.h"

using namespace ROOT::Minuit2;

struct fitParams
{
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > > Mags;
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > > Phis;
  std::map<int, double> Masses;
  std::map<int, double> Widths;
  std::map<int, double> gFactors;
  std::map<int, double> otherParams;
};


class FitParamsBase{

public:
  FitParamsBase(); 
  FitParamsBase(fitParams& theStartparams, fitParams& theErrorparams); 

  virtual ~FitParamsBase();

  std::vector< boost::shared_ptr<const JPCLS> >  jpclsVec(int index);
  virtual const std::string ampName(int index)=0;
  virtual const std::string massName(int index)=0;
  virtual const std::string widthName(int index)=0;
  virtual const std::string gFactorName(int index)=0;
  virtual const std::string otherName(int index)=0;

  virtual int ampIdxMin()=0;
  virtual int ampIdxMax()=0;
  virtual int massIdxMin()=0;
  virtual int massIdxMax()=0;
  virtual int gFactorIdxMin()=0;
  virtual int gFactorIdxMax()=0;
  virtual int otherIdxMin()=0;
  virtual int otherIdxMax()=0;

  virtual void setMnUsrParams(MnUserParameters& upar);
  virtual fitParams getFitParamVal(const std::vector<double>& par);

  int setAmpParams(std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >& ampMap, const std::vector<double>& par, int counter, int index);

  int setSingleParams(std::map<int, double>& theMap, const std::vector<double>& par, int counter, int index);

  virtual void printParams(fitParams& theParams);
  virtual void dumpParams(std::ostream& os, fitParams& theVals,  fitParams& theErrs);

protected:
  std::map < int, std::vector< boost::shared_ptr<const JPCLS> > > _jpclsMap;
  fitParams _startParams;
  fitParams _errorParams;

  std::vector<int> _ampIdx;
  std::vector<int> _massIdx;
  std::vector<int> _widthIdx;
  std::vector<int> _gFactorIdx;
  std::vector<int> _otherIdx;

private:

};


#endif
