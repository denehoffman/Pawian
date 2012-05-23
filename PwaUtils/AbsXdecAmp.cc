#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/AbsXdecAmp.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

AbsXdecAmp::AbsXdecAmp(const std::string& name, const std::vector<std::string>& hypVec) :
  _name(name)
  ,_hypVec(hypVec)
{
}

AbsXdecAmp::~AbsXdecAmp()
{
}




