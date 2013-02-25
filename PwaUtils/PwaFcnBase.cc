//#include <getopt.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <stdio.h>

#include "Minuit2/MnUserParameters.h"

#include "PwaUtils/PwaFcnBase.hh"
#include "PwaUtils/AbsLh.hh"
// #include "PwaUtils/FitParamsBase.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

PwaFcnBase::PwaFcnBase(boost::shared_ptr<AbsLh> absLh, boost::shared_ptr<FitParamsBase> fitParamsBase, std::string suffix) :
  _absLhPtr(absLh)
  , _fitParamsBasePtr(fitParamsBase)
  , _fcnCounter(0)
  , _currentResFileName("currentResult"+suffix+".dat")
{
   if (0==_absLhPtr) { Alert << "AbsLh* _absLhPtr pointer is 0 !!!!" << endmsg; exit(1); }
   _absLhPtr->getDefaultParams(_defaultFitValParms, _defaultFitErrParms);
  
}

PwaFcnBase::~PwaFcnBase()
{
}

double PwaFcnBase::operator()(const std::vector<double>& par) const
{

  fitParams theFitParmValTmp=_defaultFitValParms;

  _fitParamsBasePtr->getFitParamVal(par, theFitParmValTmp);
    
  double result=_absLhPtr->calcLogLh(theFitParmValTmp);

  _fcnCounter++;
  if (  _fcnCounter%100 == 0) {
    _fitParamsBasePtr->printParams(theFitParmValTmp);
  }
  
  if (  _fcnCounter%200 == 0) {
    std::ofstream theStream (_currentResFileName.c_str());
    _fitParamsBasePtr->dumpParams(theStream, theFitParmValTmp, (fitParams&)_defaultFitErrParms);
  }

  return result;
}

double PwaFcnBase::Up() const 
{
return .5;
}

