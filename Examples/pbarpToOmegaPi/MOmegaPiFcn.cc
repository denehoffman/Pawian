//#include <getopt.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <stdio.h>
#include "Minuit2/MnPrint.h"

#include "Examples/pbarpToOmegaPi/MOmegaPiFcn.hh"
#include "Examples/pbarpToOmegaPi/AbsOmegaPiLh.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

MOmegaPiFcn::MOmegaPiFcn(boost::shared_ptr<AbsOmegaPiLh> absOmegaPiLh) :
  _omegaPiLhPtr(absOmegaPiLh)
{
  if (0==_omegaPiLhPtr) { Alert << "OmegaPiLh pointer is 0 !!!!" << endmsg; exit(1); }
  
}

MOmegaPiFcn::~MOmegaPiFcn()
{
}

double MOmegaPiFcn::operator()(const std::vector<double>& par) const
{
  OmegaPiData::fitParamVal theFitParmValTmp;
  _omegaPiLhPtr->getFitParamVal(theFitParmValTmp, par);
 
  double result=_omegaPiLhPtr->calcLogLh(theFitParmValTmp);

//   //  print fit paramss
//   _omegaPiLhPtr->printFitParams(std::cout, theFitParmValTmp);

  return result;
}

double MOmegaPiFcn::Up() const 
{
return .5;
}

