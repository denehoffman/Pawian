//************************************************************************//
//                                                                        //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)                      //
//               Julian Pychy (julian@ep1.rub.de)                         //
//               - Ruhr-Universität Bochum                                //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

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

// The SDME for a linearly polarized photon beam
//
// rho_{l,l'} =
//
//        l = -1          l = +1
// 1/2 * (1               -f * Exp(-2ia))    l' = -1
// 1/2 * (-f * Exp(2ia)   1             )    l' = +1
//
// where "f" is the polarization fraction and "a" is the polarization angle
std::complex<double>
GammapBeamPolarization::rho(const Spin &lambda, const Spin &lambdaPrime,
                            double angle) const {
  if (lambda == lambdaPrime)
    return std::complex<double>(0.5, 0.);

  const std::complex<double> offDiagonal =
      -0.5 * _fraction * std::exp(std::complex<double>(0., -2. * angle));
  if (lambda == Spin(1) && lambdaPrime == Spin(-1))
    return offDiagonal;
  if (lambda == Spin(-1) && lambdaPrime == Spin(1))
    return std::conj(offDiagonal);

  Alert << "Photon helicities must be +1 or -1, got " << lambda << " and "
        << lambdaPrime << endmsg;
  exit(1);
}

// Calculates the intensity given two helicity amplitudes
// H_{l} and H_{l'}
//
// I = 2 * Re[sum_{l,l'} rho_{l,l'} * H_{l} * conj(H_{l'})]
double GammapBeamPolarization::intensity(
    const std::array<std::complex<double>, 2> &helicityAmps,
    double angle) const {
  const std::array<Spin, 2> helicities = {Spin(1), Spin(-1)};
  std::complex<double> result(0., 0.);

  for (unsigned int i = 0; i < helicities.size(); ++i) {
    for (unsigned int j = 0; j < helicities.size(); ++j) {
      result += rho(helicities.at(i), helicities.at(j), angle) *
                helicityAmps.at(i) * std::conj(helicityAmps.at(j));
    }
  }

  return 2. * result.real();
}

double GammapBeamPolarization::productionPlaneAngle(
    const Vector4<double> &beam, const Vector4<double> &recoil) const {
  const double beamMagnitude = beam.P();
  const double recoilMagnitude = recoil.P();
  if (beamMagnitude == 0. || recoilMagnitude == 0.)
    return _angle;

  const double beamX = beam.Px() / beamMagnitude;
  const double beamY = beam.Py() / beamMagnitude;
  const double beamZ = beam.Pz() / beamMagnitude;
  const double recoilX = recoil.Px() / recoilMagnitude;
  const double recoilY = recoil.Py() / recoilMagnitude;
  const double recoilZ = recoil.Pz() / recoilMagnitude;

  double normalX = beamZ * recoilY - beamY * recoilZ;
  double normalY = beamX * recoilZ - beamZ * recoilX;
  double normalZ = beamY * recoilX - beamX * recoilY;
  const double normalMagnitude =
      std::sqrt(normalX * normalX + normalY * normalY + normalZ * normalZ);
  if (normalMagnitude == 0.)
    return _angle;
  normalX /= normalMagnitude;
  normalY /= normalMagnitude;
  normalZ /= normalMagnitude;

  const double epsilonX = std::cos(_angle);
  const double epsilonY = std::sin(_angle);
  const double numerator = normalX * epsilonX + normalY * epsilonY;
  const double epsilonCrossNormalX = epsilonY * normalZ;
  const double epsilonCrossNormalY = -epsilonX * normalZ;
  const double epsilonCrossNormalZ =
      epsilonX * normalY - epsilonY * normalX;
  const double denominator = beamX * epsilonCrossNormalX +
                             beamY * epsilonCrossNormalY +
                             beamZ * epsilonCrossNormalZ;
  return std::atan2(numerator, denominator);
}
