//#include <getopt.h>
//#include <fstream>
//#include <string>

#include "Examples/EtacToapi0Fit/MEtacToapi0Fcn.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Data.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Lh.hh"
#include "ErrLogger/ErrLineLog.hh"
#include <cassert>

using namespace ROOT::Minuit2;

MEtacToapi0Fcn::MEtacToapi0Fcn(EtacToapi0Lh* etacToapi0Lh) :
  _etacToapi0Lh(etacToapi0Lh)
{
  if (0==_etacToapi0Lh) ErrMsg(fatal) << "EtacToapi0Lh pointer is 0 !!!!" << endmsg; 
  
}

MEtacToapi0Fcn::~MEtacToapi0Fcn()
{
}

double MEtacToapi0Fcn::operator()(const std::vector<double>& par) const
{
  fitParamVal theFitParmValTmp;
  assert(_etacToapi0Lh->setFitParamVal(theFitParmValTmp, par));
 
  double result=_etacToapi0Lh->calcLogLh(theFitParmValTmp);
  
  ErrMsg(debugging) << "InterMassFit= " << theFitParmValTmp.aMass 
		    << "  InterWidthFit= " << theFitParmValTmp.aWidth
		    << "  spin0= " << theFitParmValTmp.cont0spin 
		    << "  spin1= " << theFitParmValTmp.cont1spin 
		    << "  spin2= " << theFitParmValTmp.cont2spin 
		    << "    logLH= " << result << endmsg; 
  
  return result;
}

double MEtacToapi0Fcn::Up() const 
{
return .5;
}



