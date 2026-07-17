#pragma once

#include <complex>
#include <map>
#include <string>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

enum class Reflectivity { Negative = -1, Positive = 1 };

enum class ProtonHelicitySector { NonFlip = 0, Flip = 1 };

struct ReflectivityWaveKey {
  std::string producedSystem;
  int J;
  int M;
  Reflectivity reflectivity;
  ProtonHelicitySector sector;

  bool operator<(const ReflectivityWaveKey &other) const;
  bool operator==(const ReflectivityWaveKey &other) const;
};

class ReflectivityBasis {
public:
  using AmplitudeMap = std::map<ReflectivityWaveKey, std::complex<double>>;

  static int sign(Reflectivity reflectivity);
  static int minusOneToM(int M);
  static std::complex<double>
  fromPhotonHelicities(const std::complex<double> &plusAtM,
                       const std::complex<double> &minusAtNegativeM, int M,
                       Reflectivity reflectivity);
  static ProtonHelicitySector sector(const Spin &initialProtonHelicity,
                                     const Spin &finalProtonHelicity);
  static double parityFactor(Reflectivity reflectivity,
                             const Spin &initialProtonHelicity,
                             const Spin &finalProtonHelicity);

  static std::complex<double> reflectivityAmplitude(
      const AmplitudeMap &amplitudes, const std::string &producedSystem, int J,
      int M, Reflectivity reflectivity, const Spin &initialProtonHelicity,
      const Spin &finalProtonHelicity);

  static std::complex<double> photonHelicityAmplitude(
      const AmplitudeMap &amplitudes, const std::string &producedSystem, int J,
      int M, const Spin &photonHelicity, const Spin &initialProtonHelicity,
      const Spin &finalProtonHelicity);

  static std::string parameterStem(const ReflectivityWaveKey &key);
};
