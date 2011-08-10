#include "Examples/Tutorial/FitIF/MinimizerInterface/MIMinuitFcn.hh"
#include "Examples/Tutorial/FitIF/MinimizerInterface/MIData.hh"
#include "ErrLogger/ErrLogger.hh"
#include <cassert>

using namespace ROOT::Minuit2;

MIMinuitFcn::MIMinuitFcn(boost::shared_ptr<MIData> myData) :
  _myDataPtr(myData)
{
  if (0==_myDataPtr) {
    Alert << "Data pointer is 0 !!!!" << endmsg; 
    exit(1);
  }
}

MIMinuitFcn::~MIMinuitFcn()
{
}

double MIMinuitFcn::operator()(const std::vector<double>& par) const
{
  double result=_myDataPtr->controlParameter(par);
  DebugMsg << "current minimized value:\t"<< result << endmsg;  
  return result;
}

double MIMinuitFcn::Up() const 
{
return 1.;
}



