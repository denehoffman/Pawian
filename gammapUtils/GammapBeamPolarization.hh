#pragma once

#include <array>
#include <complex>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class GammapBeamPolarization {
public:
  GammapBeamPolarization(double fraction = 0., double angle = 0.);

  double fraction() const { return _fraction; }
  double angle() const { return _angle; }
  std::complex<double> rho(const Spin &lambda, const Spin &lambdaPrime) const;
  double
  intensity(const std::array<std::complex<double>, 2> &helicityAmps) const;

  static unsigned int helicityIndex(const Spin &lambda);

private:
  double _fraction;
  double _angle;
};
