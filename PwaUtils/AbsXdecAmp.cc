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
  ,_oldXMass(0.)
  ,_oldXWidth(0.)
  ,_oldgFactorPhiPhi(0.)
  ,_oldgFactorOmegaPhi(0.)
  ,_cacheAmps(false)
  ,_recalculate(true)
{
}

AbsXdecAmp::~AbsXdecAmp()
{
}

void AbsXdecAmp::checkRecalculation(fitParamsNew& theParamVal){
  _recalculate=true;
}


