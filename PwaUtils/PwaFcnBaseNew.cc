//#include <getopt.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <stdio.h>

#include "Minuit2/MnUserParameters.h"

#include "PwaUtils/PwaFcnBaseNew.hh"
#include "PwaUtils/AbsLhNew.hh"
// #include "PwaUtils/FitParamsBaseNew.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

PwaFcnBaseNew::PwaFcnBaseNew(boost::shared_ptr<AbsLhNew> absLh, boost::shared_ptr<FitParamsBaseNew> fitParamsBase) :
  _absLhPtr(absLh)
  , _fitParamsBasePtr(fitParamsBase)
  , _fcnCounter(0)
{
   if (0==_absLhPtr) { Alert << "AbsLh* _absLhPtr pointer is 0 !!!!" << endmsg; exit(1); }
   _absLhPtr->getDefaultParams(_defaultFitValParms, _defaultFitErrParms);
  
}

PwaFcnBaseNew::~PwaFcnBaseNew()
{
}

double PwaFcnBaseNew::operator()(const std::vector<double>& par) const
{
  _fcnCounter++;


  fitParamsNew theFitParmValTmp=_defaultFitValParms;
//   fitParamsNew theFitParmValTmp;
//   fitParamsNew theFitParmValTmp1;
//   _absLhPtr->getDefaultParams(theFitParmValTmp, theFitParmValTmp1);

  _fitParamsBasePtr->getFitParamVal(par, theFitParmValTmp);  
  
  double result=_absLhPtr->calcLogLh(theFitParmValTmp);

  DebugMsg << "logLh= " << result <<endmsg;
 
  if (  _fcnCounter%50 == 0) {  
    _fitParamsBasePtr->printParams(theFitParmValTmp);
  }
  
  if (  _fcnCounter%200 == 0) {
    std::ofstream theStream ( "currentResult.dat");
    _fitParamsBasePtr->dumpParams(theStream, theFitParmValTmp, theFitParmValTmp);
  }

  return result;
}

double PwaFcnBaseNew::Up() const 
{
return .5;
}

