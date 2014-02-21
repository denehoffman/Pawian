#pragma once

#include <cassert>
#include <vector>
#include <string>

class Particle;
class Spin;

namespace FunctionUtils {

  std::string particleListName(std::vector<Particle*>& particleVec);
  unsigned int spin3Index(Spin& spin1, Spin& spin2, Spin& spin3);

}
