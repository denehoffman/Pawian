//#include <getopt.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <stdio.h>
#include "Minuit2/MnPrint.h"

#include "Examples/pbarpToOmegaPiLS/MOmegaPiFcnLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiLhLS.hh"
#include "Examples/pbarpToOmegaPiLS/OmegaPiDataLS.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

MOmegaPiFcnLS::MOmegaPiFcnLS(boost::shared_ptr<AbsOmegaPiLhLS> absOmegaPiLh) :
  _omegaPiLhPtr(absOmegaPiLh)
{
  if (0==_omegaPiLhPtr) { Alert << "OmegaPiLh pointer is 0 !!!!" << endmsg; exit(1); }
  
}

MOmegaPiFcnLS::~MOmegaPiFcnLS()
{
}

double MOmegaPiFcnLS::operator()(const std::vector<double>& par) const
{
  OmegaPiDataLS::fitParamVal theFitParmValTmp;
  _omegaPiLhPtr->getFitParamVal(theFitParmValTmp, par);
 
  double result=_omegaPiLhPtr->calcLogLh(theFitParmValTmp);

//   //  print fit paramss
//   _omegaPiLhPtr->printFitParams(std::cout, theFitParmValTmp);

  return result;
}

double MOmegaPiFcnLS::Up() const 
{
return .5;
}

