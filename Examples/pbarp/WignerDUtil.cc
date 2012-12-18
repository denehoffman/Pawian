#include <getopt.h>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "Examples/pbarp/WignerDUtil.hh"
#include "Examples/pbarp/IsobarDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/FunctionUtils.hh"

WignerDUtil::WignerDUtil(std::vector<Particle*>& mother, std::vector<Particle*>& daughter){
  _name=FunctionUtils::particleListName(daughter)+"_"+FunctionUtils::particleListName(mother);
}

WignerDUtil::~WignerDUtil(){
}

