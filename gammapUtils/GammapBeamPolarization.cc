#include "gammapUtils/GammapBeamPolarization.hh"

#include <cmath>

#include "ErrLogger/ErrLogger.hh"

GammapBeamPolarization::GammapBeamPolarization(double fraction, double angle)
    : _fraction(fraction), _angle(angle) {}

unsigned int GammapBeamPolarization::helicityIndex(const Spin &lambda) {
  if (lambda == Spin(1))
    return 0;
  if (lambda == Spin(-1))
    return 1;
  Alert << "Photon helicity must be +1 or -1, got " << lambda << endmsg;
  exit(1);
}

std::complex<double>
GammapBeamPolarization::rho(const Spin &lambda, const Spin &lambdaPrime) const {
  if (lambda == lambdaPrime)
    return std::complex<double>(0.5, 0.);

  const std::complex<double> offDiagonal =
      -0.5 * _fraction * std::exp(std::complex<double>(0., -2. * _angle));
  if (lambda == Spin(1) && lambdaPrime == Spin(-1))
    return offDiagonal;
  if (lambda == Spin(-1) && lambdaPrime == Spin(1))
    return std::conj(offDiagonal);

  Alert << "Photon helicities must be +1 or -1, got " << lambda << " and "
        << lambdaPrime << endmsg;
  exit(1);
}

double GammapBeamPolarization::intensity(
    const std::array<std::complex<double>, 2> &helicityAmps) const {
  const std::array<Spin, 2> helicities = {Spin(1), Spin(-1)};
  std::complex<double> result(0., 0.);

  for (unsigned int i = 0; i < helicities.size(); ++i) {
    for (unsigned int j = 0; j < helicities.size(); ++j) {
      result += rho(helicities.at(i), helicities.at(j)) * helicityAmps.at(i) *
                std::conj(helicityAmps.at(j));
    }
  }

  return 2. * result.real();
}
