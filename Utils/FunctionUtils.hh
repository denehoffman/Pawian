#pragma once

#include <cassert>
#include <vector>
#include <string>
#include <sys/types.h>
#include <cstdint>

class Particle;
class Spin;

typedef uint8_t Id1StringType;
typedef uint16_t Id2StringType;
typedef uint16_t Id3StringType;

namespace FunctionUtils {

  std::string particleListName(std::vector<Particle*>& particleVec);
  Id1StringType spin1Index(const Spin& spin1);
  Id2StringType spin2Index(const Spin& spin1, const Spin& spin2);
  Id3StringType spin3Index(const Spin& spin1, const Spin& spin2, const Spin& spin3);
  short spin1IdIndex(unsigned short, const Spin& spin1);
}
