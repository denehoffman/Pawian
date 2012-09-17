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
{
}

AbsXdecAmp::~AbsXdecAmp()
{
}




