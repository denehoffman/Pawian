#pragma once

#include <math.h>

namespace PawianConstants {
  //mass of pi+- (GeV/c^2)                      
  const double mPi              = 0.13957018; 
  //mass of pi0 (GeV/c^2)
  const double mPi0             = 0.1349766;

  //square of pi+- mass 
  const double mPiSq            = mPi*mPi;
  //square of pi0 mass 
  const double mPi0Sq           = mPi0*mPi0;

  //square of pi 
  const double pi           = M_PI;

  //rad to deg 
  const double radToDeg = 180.0/M_PI;
}
