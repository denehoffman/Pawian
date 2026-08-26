#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ReflectivityBasisTest

#include <boost/test/unit_test.hpp>

#include "gammapUtils/GammapBeamPolarization.hh"
#include "gammapUtils/ReflectivityBasis.hh"

namespace {
void checkComplex(const std::complex<double> &actual,
                  const std::complex<double> &expected) {
  BOOST_CHECK_SMALL(std::abs(actual - expected), 1.e-12);
}
} // namespace

BOOST_AUTO_TEST_CASE(reflectivity_transform_round_trips) {
  const std::complex<double> positive(2.0, -0.5);
  const std::complex<double> negative(-0.25, 1.5);
  const int M = 1;

  const std::complex<double> photonPlus = negative + positive;
  const std::complex<double> photonMinusAtNegativeM =
      double(ReflectivityBasis::minusOneToM(M)) * (negative - positive);

  checkComplex(
      ReflectivityBasis::fromPhotonHelicities(
          photonPlus, photonMinusAtNegativeM, M, Reflectivity::Positive),
      positive);
  checkComplex(
      ReflectivityBasis::fromPhotonHelicities(
          photonPlus, photonMinusAtNegativeM, M, Reflectivity::Negative),
      negative);
}

BOOST_AUTO_TEST_CASE(photon_helicity_inverse_uses_negative_m) {
  ReflectivityBasis::AmplitudeMap amplitudes;
  amplitudes[{"X", 2, -1, Reflectivity::Negative,
              ProtonHelicitySector::NonFlip}] = {3., 1.};
  amplitudes[{"X", 2, -1, Reflectivity::Positive,
              ProtonHelicitySector::NonFlip}] = {1., -2.};
  amplitudes[{"X", 2, 1, Reflectivity::Negative,
              ProtonHelicitySector::NonFlip}] = {-1., 4.};
  amplitudes[{"X", 2, 1, Reflectivity::Positive,
              ProtonHelicitySector::NonFlip}] = {2., 3.};

  checkComplex(ReflectivityBasis::photonHelicityAmplitude(
                   amplitudes, "X", 2, 1, Spin(1), Spin(0.5), Spin(0.5)),
               {1., 7.});
  checkComplex(ReflectivityBasis::photonHelicityAmplitude(
                   amplitudes, "X", 2, 1, Spin(-1), Spin(0.5), Spin(0.5)),
               {-2., -3.});
}

BOOST_AUTO_TEST_CASE(parity_reduces_proton_helicity_pairs) {
  BOOST_CHECK_EQUAL(ReflectivityBasis::parityFactor(Reflectivity::Positive,
                                                    Spin(-0.5), Spin(-0.5)),
                    1.);
  BOOST_CHECK_EQUAL(ReflectivityBasis::parityFactor(Reflectivity::Negative,
                                                    Spin(-0.5), Spin(-0.5)),
                    -1.);
  BOOST_CHECK_EQUAL(ReflectivityBasis::parityFactor(Reflectivity::Positive,
                                                    Spin(-0.5), Spin(0.5)),
                    -1.);
  BOOST_CHECK_EQUAL(ReflectivityBasis::parityFactor(Reflectivity::Negative,
                                                    Spin(-0.5), Spin(0.5)),
                    1.);
}

BOOST_AUTO_TEST_CASE(parameter_names_are_deterministic) {
  const ReflectivityWaveKey key{"f2(1270)", 2, -1, Reflectivity::Positive,
                                ProtonHelicitySector::Flip};
  BOOST_CHECK_EQUAL(ReflectivityBasis::parameterStem(key),
                    "Reflectivity_f2(1270)_J2_M-1_R+_K1");
}

BOOST_AUTO_TEST_CASE(opposite_reflectivities_do_not_interfere) {
  ReflectivityBasis::AmplitudeMap amplitudes;
  amplitudes[{"X", 0, 0, Reflectivity::Negative,
              ProtonHelicitySector::NonFlip}] = {0.5, 1.25};
  amplitudes[{"X", 0, 0, Reflectivity::Positive,
              ProtonHelicitySector::NonFlip}] = {2.0, -0.75};
  amplitudes[{"X", 0, 0, Reflectivity::Negative, ProtonHelicitySector::Flip}] =
      {-1.5, 0.25};
  amplitudes[{"X", 0, 0, Reflectivity::Positive, ProtonHelicitySector::Flip}] =
      {0.4, 2.1};

  ReflectivityBasis::AmplitudeMap negativeOnly;
  ReflectivityBasis::AmplitudeMap positiveOnly;
  for (const auto &wave : amplitudes) {
    if (wave.first.reflectivity == Reflectivity::Negative)
      negativeOnly.insert(wave);
    else
      positiveOnly.insert(wave);
  }

  const GammapBeamPolarization beam(0.73, 0.31);
  const std::array<Spin, 2> protonHelicities = {Spin(0.5), Spin(-0.5)};
  const auto intensity = [&](const ReflectivityBasis::AmplitudeMap &waves) {
    double result = 0.;
    for (const Spin &initial : protonHelicities) {
      for (const Spin &final : protonHelicities) {
        const std::array<std::complex<double>, 2> photonAmplitudes = {
            ReflectivityBasis::photonHelicityAmplitude(waves, "X", 0, 0,
                                                       Spin(1), initial, final),
            ReflectivityBasis::photonHelicityAmplitude(
                waves, "X", 0, 0, Spin(-1), initial, final)};
        result += beam.intensity(photonAmplitudes, beam.angle());
      }
    }
    return result;
  };

  BOOST_CHECK_SMALL(intensity(amplitudes) - intensity(negativeOnly) -
                        intensity(positiveOnly),
                    1.e-12);
}

BOOST_AUTO_TEST_CASE(polarized_s_wave_tracks_the_production_plane) {
  const double polarization = 0.6;
  const GammapBeamPolarization beamPolarization(polarization, 0.);
  const Vector4<double> beam(8., 0., 0., 8.);
  const std::array<std::complex<double>, 2> positiveReflectivityS = {
      std::complex<double>(1., 0.), std::complex<double>(-1., 0.)};

  const Vector4<double> recoilInPolarizationPlane(2., 1., 0., 1.);
  const double parallelAngle = beamPolarization.productionPlaneAngle(
      beam, recoilInPolarizationPlane);
  BOOST_CHECK_CLOSE(
      beamPolarization.intensity(positiveReflectivityS, parallelAngle),
      2. * (1. + polarization), 1.e-10);

  const Vector4<double> recoilPerpendicularToPolarization(2., 0., 1., 1.);
  const double perpendicularAngle = beamPolarization.productionPlaneAngle(
      beam, recoilPerpendicularToPolarization);
  BOOST_CHECK_CLOSE(
      beamPolarization.intensity(positiveReflectivityS, perpendicularAngle),
      2. * (1. - polarization), 1.e-10);
}
