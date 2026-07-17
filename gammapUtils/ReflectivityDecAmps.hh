#pragma once

#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "PwaUtils/AbsXdecAmp.hh"
#include "gammapUtils/ReflectivityBasis.hh"

class AbsDecay;
class AbsPawianParameters;

class ReflectivityDecAmps : public AbsXdecAmp {
public:
  ReflectivityDecAmps(std::shared_ptr<AbsDecay> decay, ChannelID channelID);
  virtual ~ReflectivityDecAmps();

  virtual std::complex<double> XdecAmp(const Spin &lamX, EvtData *data,
                                       AbsXdecAmp *grandmaAmp = nullptr);
  virtual std::complex<double> XdecPartAmp(const Spin &lamX, Spin &lamDec,
                                           short fixDaughterNr, EvtData *data,
                                           Spin &lamFs,
                                           AbsXdecAmp *grandmaAmp = nullptr);

  std::complex<double> amplitude(const Spin &photonHelicity,
                                 const Spin &initialProtonHelicity,
                                 EvtData *data);
  std::complex<double>
  photonHelicityComponent(int M, const Spin &photonHelicity,
                          const Spin &initialProtonHelicity,
                          const Spin &finalProtonHelicity) const;
  std::complex<double>
  reflectivityComponent(int M, Reflectivity reflectivity,
                        const Spin &initialProtonHelicity,
                        const Spin &finalProtonHelicity) const;

  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();
  virtual void print(std::ostream &os) const;

  const std::vector<ReflectivityWaveKey> &waveKeys() const { return _waveKeys; }
  const std::string &producedSystem() const { return _producedSystem; }
  int producedSpin() const { return _producedSpin; }

private:
  std::string magnitudeName(const ReflectivityWaveKey &key) const;
  std::string phaseName(const ReflectivityWaveKey &key) const;
  std::complex<double> productionDynamics() const;

  std::string _producedSystem;
  int _producedSpin;
  std::vector<ReflectivityWaveKey> _waveKeys;
  ReflectivityBasis::AmplitudeMap _currentAmplitudes;
};
