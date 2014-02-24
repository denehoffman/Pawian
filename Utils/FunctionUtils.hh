#pragma once

#include <cassert>
#include <vector>
#include <string>
#include <sys/types.h>

class Particle;
class Spin;

typedef uint16_t Id2StringType;
typedef uint16_t Id3StringType;

namespace FunctionUtils {

  std::string particleListName(std::vector<Particle*>& particleVec);
  Id2StringType spin2Index(Spin& spin1, Spin& spin2);
  Id3StringType spin3Index(Spin& spin1, Spin& spin2, Spin& spin3);

}
