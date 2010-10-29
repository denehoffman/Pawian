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
{
   if (0==_psi2SToKpKmPiGamLhPtr) { Alert << "AbsPsi2SToKpKmPiGamLh* _psi2SToKpKmPiGamLhPtr pointer is 0 !!!!" << endmsg; exit(1); }
  
}

MPsi2SToKpKmPiGamFcn::~MPsi2SToKpKmPiGamFcn()
{
}

double MPsi2SToKpKmPiGamFcn::operator()(const std::vector<double>& par) const
{
  Psi2SToKpKmPiGamData::fitParamVal theFitParmValTmp;
  _psi2SToKpKmPiGamLhPtr->setFitParamVal(theFitParmValTmp, par);
 
  double result=_psi2SToKpKmPiGamLhPtr->calcLogLh(theFitParmValTmp);

  DebugMsg << "logLh= " << result <<endmsg;  
  _psi2SToKpKmPiGamLhPtr->printCurrentFitResult(theFitParmValTmp);

  return result;
}

double MPsi2SToKpKmPiGamFcn::Up() const 
{
return .5;
}

