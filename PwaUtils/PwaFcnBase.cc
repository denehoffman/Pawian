//#include <getopt.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <stdio.h>

#include "Minuit2/MnUserParameters.h"

#include "PwaUtils/PwaFcnBase.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

PwaFcnBase::PwaFcnBase(boost::shared_ptr<AbsLh> absLh, boost::shared_ptr<FitParamsBase> fitParamsBase) :
  _absLhPtr(absLh)
  , _fitParamsBasePtr(fitParamsBase)
  , _fcnCounter(new unsigned int (0))
{
   if (0==_absLhPtr) { Alert << "AbsLh* _absLhPtr pointer is 0 !!!!" << endmsg; exit(1); }
  
}

PwaFcnBase::~PwaFcnBase()
{
  delete _fcnCounter;
}

double PwaFcnBase::operator()(const std::vector<double>& par) const
{
  (*_fcnCounter)++;
  fitParams theFitParmValTmp=_fitParamsBasePtr->getFitParamVal(par);
 
 
  
  double result=_absLhPtr->calcLogLh(theFitParmValTmp);

  DebugMsg << "logLh= " << result <<endmsg;
 
  if (  (*_fcnCounter)%10 == 0) {  
    _fitParamsBasePtr->printParams(theFitParmValTmp);
    for(int i=0;i<par.size();i++){
      Info << "Fit parameter " << i << " " << par[i] << endmsg;
    }
  }
  
  if (  (*_fcnCounter)%100 == 0) {
    std::ofstream theStream ( "currentResult.dat");
    _fitParamsBasePtr->dumpParams(theStream, theFitParmValTmp, theFitParmValTmp);
  }

  return result;
}

double PwaFcnBase::Up() const 
{
return .5;
}

