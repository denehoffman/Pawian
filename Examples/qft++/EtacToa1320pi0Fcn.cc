//#include <getopt.h>
//#include <fstream>
//#include <string>

#include "Examples/qft++/EtacToa1320pi0Fcn.hh"
#include "Examples/qft++/EtacToa1320pi0Data.hh"
#include "Examples/qft++/EtacToa1320pi0fit.hh"
#include "ErrLogger/ErrLineLog.hh"
#include <cassert>

using namespace ROOT::Minuit2;

EtacToa1320pi0Fcn::EtacToa1320pi0Fcn(EtacToa1320pi0fit* etacToa1320pi0fit) :
  _etacToa1320pi0fit(etacToa1320pi0fit)
{
  if (0==_etacToa1320pi0fit) ErrMsg(fatal) << "EtacToa1320pi0fit pointer is 0 !!!!" << endmsg; 
  
}

EtacToa1320pi0Fcn::~EtacToa1320pi0Fcn()
{
}

double EtacToa1320pi0Fcn::operator()(const std::vector<double>& par) const
{
  fitParamVal theFitParmValTmp;
  assert(_etacToa1320pi0fit->setFitParamVal(theFitParmValTmp, par));
 
  double result=_etacToa1320pi0fit->calcLogLh(theFitParmValTmp);
  
  ErrMsg(debugging) << "InterMassFit= " << theFitParmValTmp.a1320Mass 
		    << "  InterWidthFit= " << theFitParmValTmp.a1320Width
		    << "  spin0= " << theFitParmValTmp.cont0spin 
		    << "  spin1= " << theFitParmValTmp.cont1spin 
		    << "  spin2= " << theFitParmValTmp.cont2spin 
		    << "    logLH= " << result << endmsg; 
  
  return result;
}

double EtacToa1320pi0Fcn::Up() const 
{
return .5;
}



