//#include <getopt.h>
//#include <fstream>
//#include <string>

#include "Examples/EtacToapi0Fit/MEtacToapi0Fcn.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Data.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Lh.hh"
#include "ErrLogger/ErrLineLog.hh"
#include <cassert>

using namespace ROOT::Minuit2;

MEtacToapi0Fcn::MEtacToapi0Fcn(boost::shared_ptr<EtacToapi0Lh> etacToapi0Lh) :
  _etacToapi0LhPtr(etacToapi0Lh)
{
  if (0==_etacToapi0LhPtr) ErrMsg(fatal) << "EtacToapi0Lh pointer is 0 !!!!" << endmsg; 
  
}

MEtacToapi0Fcn::~MEtacToapi0Fcn()
{
}

double MEtacToapi0Fcn::operator()(const std::vector<double>& par) const
{
  fitParamVal theFitParmValTmp;
  assert(_etacToapi0LhPtr->setFitParamVal(theFitParmValTmp, par));
 
  double result=_etacToapi0LhPtr->calcLogLh(theFitParmValTmp);
  
  ErrMsg(debugging) << "InterMassFit= " << theFitParmValTmp.aMass << "\n" 
		    << "  InterWidthFit= " << theFitParmValTmp.aWidth << "\n"
		    << "  spin0= " << theFitParmValTmp.cont0spin << "\n" 
		    << "  spin1= " << theFitParmValTmp.cont1spin << "\n" 
		    << "  spin2= " << theFitParmValTmp.cont2spin << "\n" 
		    << "    logLH= " << result << endmsg; 
  
  return result;
}

double MEtacToapi0Fcn::Up() const 
{
return .5;
}



