//#include <getopt.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Minuit2/MnUserParameters.h"

#include "PwaUtils/PwaFcnBaseNew.hh"
#include "PwaUtils/AbsLhNew.hh"
// #include "PwaUtils/FitParamsBaseNew.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

PwaFcnBaseNew::PwaFcnBaseNew(boost::shared_ptr<AbsLhNew> absLh, boost::shared_ptr<FitParamsBaseNew> fitParamsBase, std::string suffix) :
  _absLhPtr(absLh)
  , _fitParamsBasePtr(fitParamsBase)
  , _fcnCounter(0)
  , _currentResFileName("currentResult"+suffix+".dat")
{
   if (0==_absLhPtr) { Alert << "AbsLh* _absLhPtr pointer is 0 !!!!" << endmsg; exit(1); }
   _absLhPtr->getDefaultParams(_defaultFitValParms, _defaultFitErrParms);
  
}

PwaFcnBaseNew::~PwaFcnBaseNew()
{
}

double PwaFcnBaseNew::operator()(const std::vector<double>& par) const
{
  

  fitParamsNew theFitParmValTmp=_defaultFitValParms;

#ifdef _OPENMP
#pragma omp critical
  {
 #endif

     _fitParamsBasePtr->getFitParamVal(par, theFitParmValTmp);
    
#ifdef _OPENMP  
  }
#endif


  double result=_absLhPtr->calcLogLh(theFitParmValTmp);


#ifdef _OPENMP
#pragma omp critical
  {
#endif
    _fcnCounter++;
    DebugMsg << "logLh= " << result <<endmsg;
    
    if (  _fcnCounter%50 == 0) {  
      _fitParamsBasePtr->printParams(theFitParmValTmp);
    }
    
    if (  _fcnCounter%200 == 0) {
      std::ofstream theStream (_currentResFileName.c_str());
      _fitParamsBasePtr->dumpParams(theStream, theFitParmValTmp, theFitParmValTmp);
    }

#ifdef _OPENMP  
  }
#endif

  return result;
}

double PwaFcnBaseNew::Up() const 
{
return .5;
}

