#pragma once

#include <complex>

namespace SphericalHarmonic{
  double LegendrePolynomial(unsigned int l, int m, double cosT);
  std::complex<double> SphericalHarmonicFunc(int l, int m, double theta, double phi);
}
