#include "gammapUtils/ReflectivityDecAmps.hh"

#include <cmath>

#include "ErrLogger/ErrLogger.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "Particle/Particle.hh"
#include "PwaUtils/AbsDecay.hh"

ReflectivityDecAmps::ReflectivityDecAmps(std::shared_ptr<AbsDecay> decay,
                                         ChannelID channelID)
    : AbsXdecAmp(decay, channelID), _producedSystem(""), _producedSpin(0) {
  AbsXdecAmp::initialize();

  if (_daughter1IsStable) {
    Alert << "Reflectivity production requires daughter 1 to be the produced "
             "decaying system"
          << endmsg;
    exit(1);
  }
  if (!_daughter2IsStable || _Jdaughter2 != Spin(0.5)) {
    Alert << "Reflectivity production requires a stable spin-1/2 recoil as "
             "daughter 2"
          << endmsg;
    exit(1);
  }
  if (int(2. * _Jdaughter1) % 2 != 0) {
    Alert << "The photoproduction reflectivity basis requires integer spin "
             "for the produced system, got J="
          << _Jdaughter1 << endmsg;
    exit(1);
  }
  if (_absDyn->isLdependent()) {
    Alert << "L-dependent root production dynamics are not supported by the "
             "reflectivity formalism"
          << endmsg;
    exit(1);
  }

  _producedSystem = _daughter1Name;
  _producedSpin = int(_Jdaughter1);
  fillParamNameList();
}

ReflectivityDecAmps::~ReflectivityDecAmps() {}

std::string
ReflectivityDecAmps::magnitudeName(const ReflectivityWaveKey &key) const {
  return ReflectivityBasis::parameterStem(key) + "Mag";
}

std::string
ReflectivityDecAmps::phaseName(const ReflectivityWaveKey &key) const {
  return ReflectivityBasis::parameterStem(key) + "Phi";
}

void ReflectivityDecAmps::fillParamNameList() {
  _paramNameList.clear();
  _waveKeys.clear();
  const Reflectivity reflectivities[] = {Reflectivity::Negative,
                                         Reflectivity::Positive};
  const ProtonHelicitySector sectors[] = {ProtonHelicitySector::NonFlip,
                                          ProtonHelicitySector::Flip};

  for (int M = -_producedSpin; M <= _producedSpin; ++M) {
    for (Reflectivity reflectivity : reflectivities) {
      for (ProtonHelicitySector sector : sectors) {
        const ReflectivityWaveKey key{_producedSystem, _producedSpin, M,
                                      reflectivity, sector};
        _waveKeys.push_back(key);
        _paramNameList.push_back(magnitudeName(key));
        _paramNameList.push_back(phaseName(key));
      }
    }
  }
}

// Fills the default params with a magnitude of 1/sqrt(N) where N is the number
// of waves
void ReflectivityDecAmps::fillDefaultParams(
    std::shared_ptr<AbsPawianParameters> fitPar) {
  const double defaultMagnitude =
      _waveKeys.empty() ? 1. : 1. / std::sqrt(double(_waveKeys.size()));
  for (const ReflectivityWaveKey &key : _waveKeys) {
    fitPar->Add(magnitudeName(key), defaultMagnitude, defaultMagnitude / 2.);
    fitPar->Add(phaseName(key), 0., 0.2);
  }

  _absDyn->fillDefaultParams(fitPar);
  if (!_daughter1IsStable)
    _decAmpDaughter1->fillDefaultParams(fitPar);
  if (!_daughter2IsStable)
    _decAmpDaughter2->fillDefaultParams(fitPar);
}

void ReflectivityDecAmps::updateFitParams(
    std::shared_ptr<AbsPawianParameters> fitPar) {
  for (const ReflectivityWaveKey &key : _waveKeys) {
    const double magnitude = std::abs(fitPar->Value(magnitudeName(key)));
    const double phase = fitPar->Value(phaseName(key));
    _currentAmplitudes[key] =
        _preFactor * _isospinCG * std::polar(magnitude, phase);
  }

  _absDyn->updateFitParams(fitPar);
  if (!_daughter1IsStable)
    _decAmpDaughter1->updateFitParams(fitPar);
  if (!_daughter2IsStable)
    _decAmpDaughter2->updateFitParams(fitPar);
}

// This is the amplitude A([R],L,M;k), see ReflectivityBasis.cc
std::complex<double> ReflectivityDecAmps::reflectivityComponent(
    int M, Reflectivity reflectivity, const Spin &initialProtonHelicity,
    const Spin &finalProtonHelicity) const {
  return ReflectivityBasis::reflectivityAmplitude(
      _currentAmplitudes, _producedSystem, _producedSpin, M, reflectivity,
      initialProtonHelicity, finalProtonHelicity);
}

// This is the amplitude T(J,M,lda_gamma;lda_p,lda_p'), see ReflectivityBasis.cc
std::complex<double> ReflectivityDecAmps::photonHelicityComponent(
    int M, const Spin &photonHelicity, const Spin &initialProtonHelicity,
    const Spin &finalProtonHelicity) const {
  return ReflectivityBasis::photonHelicityAmplitude(
      _currentAmplitudes, _producedSystem, _producedSpin, M, photonHelicity,
      initialProtonHelicity, finalProtonHelicity);
}

std::complex<double> ReflectivityDecAmps::productionDynamics() const {
  return _cachedDynIdMap.at(_absDyn->grandMaId(nullptr));
}

//
std::complex<double>
ReflectivityDecAmps::amplitude(const Spin &photonHelicity,
                               const Spin &initialProtonHelicity,
                               EvtData *data) {
  const Spin finalProtonHelicity = _lam2MinProj;
  std::complex<double> result(0., 0.);
  for (int M = -_producedSpin; M <= _producedSpin; ++M) {
    result += photonHelicityComponent(M, photonHelicity, initialProtonHelicity,
                                      finalProtonHelicity) *
              daughterAmp(Spin(M), finalProtonHelicity, data);
  }
  return result * productionDynamics();
}

std::complex<double> ReflectivityDecAmps::XdecAmp(const Spin &, EvtData *,
                                                  AbsXdecAmp *) {
  Alert << "ReflectivityDecAmps requires separate photon and proton "
           "helicities; use amplitude()"
        << endmsg;
  exit(1);
}

std::complex<double> ReflectivityDecAmps::XdecPartAmp(const Spin &, Spin &,
                                                      short, EvtData *, Spin &,
                                                      AbsXdecAmp *) {
  Alert << "ReflectivityDecAmps does not support the XdecPartAmp "
           "interface"
        << endmsg;
  exit(1);
}

void ReflectivityDecAmps::print(std::ostream &os) const {
  os << "Reflectivity amplitudes for " << _producedSystem
     << " (J=" << _producedSpin << ", " << _waveKeys.size() << " waves)";
}
