#include "Examples/Tutorial/MinuitFit/MinuitFitFcn.hh"
#include "Examples/Tutorial/MinuitFit/MinuitFit.hh"
#include "ErrLogger/ErrLogger.hh"
#include <cassert>

using namespace ROOT::Minuit2;

MinuitFitFcn::MinuitFitFcn(boost::shared_ptr<MinuitFit> minuitFit) :
  _minFitPtr(minuitFit)
{
  if (0==_minFitPtr) {
    Alert << "MinuitFit pointer is 0 !!!!" << endmsg; 
    exit(1);
  }
}

MinuitFitFcn::~MinuitFitFcn()
{
}

double MinuitFitFcn::operator()(const std::vector<double>& par) const
{
  double result=_minFitPtr->calcChiSqr(par);
  DebugMsg << "current chi^2:\t"<< result << endmsg;  
  return result;
}

double MinuitFitFcn::Up() const 
{
return 1.;
}



