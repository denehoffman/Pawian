#pragma once

#include <math.h>
#include <complex>

namespace PawianConstants {

  //mass of electron (GeV/c^2)                      
  const double mElectron           = 0.0005109989;
  //mass of proton (GeV/c^2)                      
  const double mProton              = 0.938272081; 
  //mass of pi+- (GeV/c^2)                      
  const double mPi              = 0.13957018; 
  //mass of pi0 (GeV/c^2)
  const double mPi0             = 0.1349766;

  //square of pi+- mass 
  const double mPiSq            = mPi*mPi;
  //square of pi0 mass 
  const double mPi0Sq           = mPi0*mPi0;

  //constant pi 
  const double pi           = M_PI;

  //constant 2*pi 
  const double twopi           = 2.*M_PI;

  //rad to deg 
  const double radToDeg = 180.0/M_PI;

  //deg to rad 
  const double degToRad = M_PI/180.0;

  //complex number i
  const std::complex<double> i(0.,1.); 
}
