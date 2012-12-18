#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/AbsXdecAmp.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

AbsXdecAmp::AbsXdecAmp(const std::string& name, const std::vector<std::string>& hypVec, Spin spinX, int parity) :
  _name(name)
  ,_hypVec(hypVec)
  ,_J_X(spinX)
  ,_parity(parity)
  ,_currentXMass(0.)
  ,_currentXWidth(0.)
  ,_currentgFactorPhiPhi(0.)
  ,_currentgFactorOmegaPhi(0.)
  ,_cacheAmps(false)
  ,_recalculate(true)
{
}

AbsXdecAmp::AbsXdecAmp(const std::string& name) :
 _name(name)
{
}

AbsXdecAmp::~AbsXdecAmp()
{
}

bool AbsXdecAmp::checkRecalculation(fitParamsNew& theParamVal){
  _recalculate=true;
  return _recalculate;
}


