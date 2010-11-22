//#include <getopt.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <stdio.h>

#include "Minuit2/MnUserParameters.h"

#include "Examples/Psi2STo2K2PiGam/MPsi2STo2K2PiGamFcn.hh"
#include "Examples/Psi2STo2K2PiGam/AbsPsi2STo2K2PiGamLh.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

MPsi2STo2K2PiGamFcn::MPsi2STo2K2PiGamFcn(boost::shared_ptr<AbsPsi2STo2K2PiGamLh> psi2STo2K2PiGamLh) :
  _psi2STo2K2PiGamLhPtr(psi2STo2K2PiGamLh),
  _fcnCounter(new unsigned int (0))
{
   if (0==_psi2STo2K2PiGamLhPtr) { Alert << "AbsPsi2STo2K2PiGamLh* _psi2STo2K2PiGamLhPtr pointer is 0 !!!!" << endmsg; exit(1); }
//    (*_fcnCounter)=0;  
}

MPsi2STo2K2PiGamFcn::~MPsi2STo2K2PiGamFcn()
{
  delete _fcnCounter;
}

double MPsi2STo2K2PiGamFcn::operator()(const std::vector<double>& par) const
{
  (*_fcnCounter)++;
  Psi2STo2K2PiGamData::fitParamVal theFitParmValTmp;
  _psi2STo2K2PiGamLhPtr->setFitParamVal(theFitParmValTmp, par);
 
  double result=_psi2STo2K2PiGamLhPtr->calcLogLh(theFitParmValTmp);

  DebugMsg << "logLh= " << result <<endmsg;  
  _psi2STo2K2PiGamLhPtr->printCurrentFitResult(theFitParmValTmp);

  if (  (*_fcnCounter)%100 == 0) {
    std::ofstream theStream ( "currentResult.dat");
    std::string theSuffix="Val"; 
    _psi2STo2K2PiGamLhPtr->dumpCurrentResult(theStream, theFitParmValTmp, theSuffix);
  }
  return result;
}

double MPsi2STo2K2PiGamFcn::Up() const 
{
return .5;
}

