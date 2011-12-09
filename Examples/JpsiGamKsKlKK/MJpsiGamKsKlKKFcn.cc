//#include <getopt.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <stdio.h>

#include "Minuit2/MnUserParameters.h"

#include "Examples/JpsiGamKsKlKK/MJpsiGamKsKlKKFcn.hh"
#include "Examples/JpsiGamKsKlKK/AbsJpsiGamKsKlKKLh.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

MJpsiGamKsKlKKFcn::MJpsiGamKsKlKKFcn(boost::shared_ptr<AbsJpsiGamKsKlKKLh> psi2SToKpKmPiGamLh) :
  _JpsiGamKsKlKKLhPtr(psi2SToKpKmPiGamLh)
  ,_fcnCounter(new unsigned int (0))
{
   if (0==_JpsiGamKsKlKKLhPtr) { Alert << "AbsJpsiGamKsKlKKLh* _JpsiGamKsKlKKLhPtr pointer is 0 !!!!" ;  // << endmsg; 
     exit(1); 
   }
  
}

MJpsiGamKsKlKKFcn::~MJpsiGamKsKlKKFcn()
{
  delete _fcnCounter;
}

double MJpsiGamKsKlKKFcn::operator()(const std::vector<double>& par) const
{
  (*_fcnCounter)++;
  paramGamKsKlKK theFitParmValTmp;
  _JpsiGamKsKlKKLhPtr->setFitParamVal(theFitParmValTmp, par);
 
  double result=_JpsiGamKsKlKKLhPtr->calcLogLh(theFitParmValTmp);

  Info << "logLh= " << result ;  // << endmsg;

  if (  (*_fcnCounter)%20 == 0) {  
  _JpsiGamKsKlKKLhPtr->printCurrentFitResult(theFitParmValTmp);
  }

  if (  (*_fcnCounter)%100 == 0) {
    std::ofstream theStream ( "currentResult.dat");
    std::string theSuffix="Val"; 
    _JpsiGamKsKlKKLhPtr->dumpCurrentResult(theStream, theFitParmValTmp, theSuffix);
  }

  return result;
}

double MJpsiGamKsKlKKFcn::Up() const 
{
return .5;
}

