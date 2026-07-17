#include "gammapUtils/ReflectivityBasis.hh"

#include <sstream>
#include <tuple>

#include "ErrLogger/ErrLogger.hh"

// https://arxiv.org/pdf/1906.04841
// We define helicity amplitudes as T(J, M, lda_gamma; lda_p, lda_p')
//
// Reflectivity amplitudes are then defined as
//
// T([R], J, M; lda_p, lda_p')
//     = (1/2)[T(J,M,+;...) - R(-1)^M T(J,-M,-;...)]
// (equation D1)
//
// There is an additional parity relation:
//
// T([R], J, M, lda_gamma; -lda_p, -lda_p')
//     = R(-1)^{lda_p - lda_p'} T([R], J, M; lda_p, lda_p')
//
// This means the only independent amplitudes are those where the spin flips
// (lda_p, lda_p') = (+-) or (-+) and those where it doesn't (++) or (--)
//
// Then we can define a simpler set of amplitudes to actually use:
//
// A([R],J,M;0) = T([R],J,M;++) = R T([R],J,M;--)
// A([R],J,M;1) = T([R],J,M;+-) = -R T([R],J,M;-+)
//
// (see parityFactor below)

bool ReflectivityWaveKey::operator<(const ReflectivityWaveKey &other) const {
  return std::tie(producedSystem, J, M, reflectivity, sector) <
         std::tie(other.producedSystem, other.J, other.M, other.reflectivity,
                  other.sector);
}

bool ReflectivityWaveKey::operator==(const ReflectivityWaveKey &other) const {
  return producedSystem == other.producedSystem && J == other.J &&
         M == other.M && reflectivity == other.reflectivity &&
         sector == other.sector;
}

int ReflectivityBasis::sign(Reflectivity reflectivity) {
  return static_cast<int>(reflectivity);
}

// Helper for (-1)^M
int ReflectivityBasis::minusOneToM(int M) { return (M % 2 == 0) ? 1 : -1; }

// Equation D1
std::complex<double> ReflectivityBasis::fromPhotonHelicities(
    const std::complex<double> &plusAtM,
    const std::complex<double> &minusAtNegativeM, int M,
    Reflectivity reflectivity) {
  return 0.5 * (plusAtM -
                double(sign(reflectivity) * minusOneToM(M)) * minusAtNegativeM);
}

// Validation for proton spin sectors
ProtonHelicitySector
ReflectivityBasis::sector(const Spin &initialProtonHelicity,
                          const Spin &finalProtonHelicity) {
  if (std::abs(initialProtonHelicity) != Spin(0.5) ||
      std::abs(finalProtonHelicity) != Spin(0.5)) {
    Alert << "Proton helicities must be +/-1/2, got " << initialProtonHelicity
          << " and " << finalProtonHelicity << endmsg;
    exit(1);
  }
  return initialProtonHelicity == finalProtonHelicity
             ? ProtonHelicitySector::NonFlip
             : ProtonHelicitySector::Flip;
}

// R | lda_p | lda_p' | factor
// + | +1/2  | +1/2   | +1
// + | +1/2  | -1/2   | +1
// - | +1/2  | +1/2   | +1
// - | +1/2  | -1/2   | +1
// + | -1/2  | +1/2   | -1
// + | -1/2  | -1/2   | +1
// - | -1/2  | +1/2   | +1
// - | -1/2  | -1/2   | -1
//
// This seems weird, but it makes sense if we write:
//
// (lda_p, lda_p') | Amplitude
// (++)            | A([R],J,M;0) = T([R],J,M;++)
// (+-)            | A([R],J,M;1) = T([R],J,M;+-)
// (--)            | R * A([R],J,M;0) = T([R],J,M;--)
// (-+)            | -R * A([R],J,M;1) = T([R],J,M;-+)
//
// All lda_p = + terms don't get any parity or reflectivity factor since we
// define the spin-flip basis as (++), (+-), but to convert (--) and (-+), we
// need additional factors to preserve the parity relations.
double ReflectivityBasis::parityFactor(Reflectivity reflectivity,
                                       const Spin &initialProtonHelicity,
                                       const Spin &finalProtonHelicity) {
  sector(initialProtonHelicity, finalProtonHelicity);
  if (initialProtonHelicity == Spin(0.5))
    return 1.;

  const double helicityPairFactor =
      initialProtonHelicity == finalProtonHelicity ? 1. : -1.;
  return sign(reflectivity) * helicityPairFactor;
}

// parityFactor * A([R],J,M;lda_p,lda_p') => T([R],J,M,k;lda_p,lda_p')
std::complex<double> ReflectivityBasis::reflectivityAmplitude(
    const AmplitudeMap &amplitudes, const std::string &producedSystem, int J,
    int M, Reflectivity reflectivity, const Spin &initialProtonHelicity,
    const Spin &finalProtonHelicity) {
  const ReflectivityWaveKey key{
      producedSystem, J, M, reflectivity,
      sector(initialProtonHelicity, finalProtonHelicity)};
  const auto it = amplitudes.find(key);
  if (it == amplitudes.end())
    return std::complex<double>(0., 0.);
  return parityFactor(reflectivity, initialProtonHelicity,
                      finalProtonHelicity) *
         it->second;
}

// This converts the reflectivity amplitudes to helicity amplitudes
// Equation D2:
// T(J,M,+1;...) = [ T([-],J,M,...) + T([+],J,M,...) ]
// T(J,M,-1;...) = (-1)^M [ T([-],J,-M;...) - T([+],J,-M;...) ]
std::complex<double> ReflectivityBasis::photonHelicityAmplitude(
    const AmplitudeMap &amplitudes, const std::string &producedSystem, int J,
    int M, const Spin &photonHelicity, const Spin &initialProtonHelicity,
    const Spin &finalProtonHelicity) {
  if (photonHelicity == Spin(1)) {
    return reflectivityAmplitude(amplitudes, producedSystem, J, M,
                                 Reflectivity::Negative, initialProtonHelicity,
                                 finalProtonHelicity) +
           reflectivityAmplitude(amplitudes, producedSystem, J, M,
                                 Reflectivity::Positive, initialProtonHelicity,
                                 finalProtonHelicity);
  }
  if (photonHelicity == Spin(-1)) {
    return static_cast<double>(minusOneToM(M)) *
           (reflectivityAmplitude(amplitudes, producedSystem, J, -M,
                                  Reflectivity::Negative, initialProtonHelicity,
                                  finalProtonHelicity) -
            reflectivityAmplitude(amplitudes, producedSystem, J, -M,
                                  Reflectivity::Positive, initialProtonHelicity,
                                  finalProtonHelicity));
  }

  Alert << "Photon helicity must be +/-1, got " << photonHelicity << endmsg;
  exit(1);
}

std::string ReflectivityBasis::parameterStem(const ReflectivityWaveKey &key) {
  std::ostringstream result;
  result << "Reflectivity_" << key.producedSystem << "_J" << key.J << "_M"
         << key.M << "_R"
         << (key.reflectivity == Reflectivity::Positive ? "+" : "-") << "_K"
         << (key.sector == ProtonHelicitySector::NonFlip ? 0 : 1);
  return result.str();
}
