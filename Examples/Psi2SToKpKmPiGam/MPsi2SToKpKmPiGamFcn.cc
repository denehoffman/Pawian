//#include <getopt.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <stdio.h>

#include "Minuit2/MnUserParameters.h"

#include "Examples/Psi2SToKpKmPiGam/MPsi2SToKpKmPiGamFcn.hh"
#include "Examples/Psi2SToKpKmPiGam/AbsPsi2SToKpKmPiGamLh.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

MPsi2SToKpKmPiGamFcn::MPsi2SToKpKmPiGamFcn(boost::shared_ptr<AbsPsi2SToKpKmPiGamLh> psi2SToKpKmPiGamLh) :
  _psi2SToKpKmPiGamLhPtr(psi2SToKpKmPiGamLh)
  ,_fcnCounter(new unsigned int (0))
{
   if (0==_psi2SToKpKmPiGamLhPtr) { Alert << "AbsPsi2SToKpKmPiGamLh* _psi2SToKpKmPiGamLhPtr pointer is 0 !!!!" << endmsg; exit(1); }
  
}

MPsi2SToKpKmPiGamFcn::~MPsi2SToKpKmPiGamFcn()
{
  delete _fcnCounter;
}

double MPsi2SToKpKmPiGamFcn::operator()(const std::vector<double>& par) const
{
  (*_fcnCounter)++;
  paramKpKmPiGam theFitParmValTmp;
  _psi2SToKpKmPiGamLhPtr->setFitParamVal(theFitParmValTmp, par);
 
  double result=_psi2SToKpKmPiGamLhPtr->calcLogLh(theFitParmValTmp);

  DebugMsg << "logLh= " << result <<endmsg;

  if (  (*_fcnCounter)%20 == 0) {  
  _psi2SToKpKmPiGamLhPtr->printCurrentFitResult(theFitParmValTmp);
  }

  if (  (*_fcnCounter)%100 == 0) {
    std::ofstream theStream ( "currentResult.dat");
    std::string theSuffix="Val"; 
    _psi2SToKpKmPiGamLhPtr->dumpCurrentResult(theStream, theFitParmValTmp, theSuffix);
  }

  return result;
}

double MPsi2SToKpKmPiGamFcn::Up() const 
{
return .5;
}

