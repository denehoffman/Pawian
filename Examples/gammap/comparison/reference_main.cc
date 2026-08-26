#include "ROOTDataReader.h"
#include "Zlm.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/NormIntInterface.h"

#include "TLorentzVector.h"

#include <cctype>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>

#ifndef AMPTOOLS_REFERENCE_REVISION
#define AMPTOOLS_REFERENCE_REVISION "unknown"
#endif

#ifndef PAWIAN_REVISION
#define PAWIAN_REVISION "unknown"
#endif

namespace {

constexpr const char *kReaction = "reference";
constexpr const char *kAmplitude = "reference::PositiveReal::S0+";

void registerComponents() {
  AmpToolsInterface::registerDataReader(ROOTDataReader());
  AmpToolsInterface::registerAmplitude(Zlm());
}

bool hasExplicitWeights(DataReader *reader) {
  auto *rootReader = dynamic_cast<ROOTDataReader *>(reader);
  return rootReader != nullptr && rootReader->hasWeight();
}

bool isSha256(const std::string &digest) {
  if (digest.size() != 64)
    return false;
  for (const char character : digest) {
    if (!std::isxdigit(static_cast<unsigned char>(character)))
      return false;
  }
  return true;
}

bool writeReaderAudit(std::ostream &output, DataReader *reader) {
  reader->resetSource();
  output << "{\"event_count\":" << reader->numEvents() << ",\"events\":[";
  for (unsigned int index = 0; index < reader->numEvents(); ++index) {
    std::unique_ptr<Kinematics> event(reader->getEvent());
    if (!event)
      return false;
    if (index != 0)
      output << ',';
    output << "{\"event_id\":" << event->eventID()
           << ",\"weight\":" << event->weight() << '}';
  }
  output << "]}";
  return true;
}

int prepareNormalization(const std::string &configPath,
                         const std::string &normalizationPath,
                         const std::string &auditPath) {
  registerComponents();
  ConfigFileParser parser(configPath);
  AmpToolsInterface interface(parser.getConfigurationInfo());
  NormIntInterface *normalization = interface.normIntInterface(kReaction);
  if (normalization == nullptr || !normalization->hasAccessToMC()) {
    std::cerr << "normalization configuration does not provide generated and "
                 "accepted MC\n";
    return 2;
  }
  if (!hasExplicitWeights(interface.genMCReader(kReaction)) ||
      !hasExplicitWeights(interface.accMCReader(kReaction))) {
    std::cerr << "normalization inputs must provide explicit Weight branches\n";
    return 2;
  }

  normalization->forceCacheUpdate();
  normalization->exportNormIntCache(normalizationPath);

  std::ofstream audit(auditPath);
  if (!audit) {
    std::cerr << "cannot open normalization audit output: " << auditPath
              << '\n';
    return 2;
  }
  audit << std::setprecision(17) << "{\"generated\":";
  if (!writeReaderAudit(audit, interface.genMCReader(kReaction))) {
    std::cerr << "generated normalization reader ended early\n";
    return 2;
  }
  audit << ",\"accepted\":";
  if (!writeReaderAudit(audit, interface.accMCReader(kReaction))) {
    std::cerr << "accepted normalization reader ended early\n";
    return 2;
  }
  audit << "}\n";
  std::cout << "prepared normalization: generated="
            << normalization->numGenEvents()
            << " accepted_weight_sum=" << normalization->numAccEvents()
            << " file=" << normalizationPath << '\n';
  return 0;
}

void writeFourMomentum(std::ostream &output, const TLorentzVector &momentum) {
  output << '[' << momentum.E() << ',' << momentum.Px() << ',' << momentum.Py()
         << ',' << momentum.Pz() << ']';
}

int evaluate(const std::string &configPath, const std::string &outputPath,
             const std::string &sampleDigest,
             const std::string &normalizationAuditPath) {
  registerComponents();
  ConfigFileParser parser(configPath);
  AmpToolsInterface interface(parser.getConfigurationInfo());
  DataReader *reader = interface.dataReader(kReaction);
  if (!hasExplicitWeights(reader)) {
    std::cerr << "data input must provide an explicit Weight branch\n";
    return 2;
  }
  interface.loadEvents(reader);
  interface.processEvents(kReaction);

  NormIntInterface *normalization = interface.normIntInterface(kReaction);
  if (normalization == nullptr || normalization->hasAccessToMC()) {
    std::cerr << "evaluation must load the explicitly prepared normalization "
                 "file\n";
    return 2;
  }

  std::ofstream output(outputPath);
  if (!output) {
    std::cerr << "cannot open reference output: " << outputPath << '\n';
    return 2;
  }
  std::ifstream auditInput(normalizationAuditPath);
  if (!auditInput) {
    std::cerr << "cannot open normalization audit: " << normalizationAuditPath
              << '\n';
    return 2;
  }
  const std::string normalizationAudit(
      (std::istreambuf_iterator<char>(auditInput)),
      std::istreambuf_iterator<char>());
  if (normalizationAudit.empty()) {
    std::cerr << "normalization audit is empty\n";
    return 2;
  }
  output << std::setprecision(17);

  const std::complex<double> norm =
      normalization->normInt(kAmplitude, kAmplitude, true);
  output << "{\n"
         << "  \"format_version\": 1,\n"
         << "  \"framework\": {\"name\": \"AmpTools\", \"revision\": \""
         << AMPTOOLS_REFERENCE_REVISION << "\"},\n"
         << "  \"pawian_revision\": \"" << PAWIAN_REVISION << "\",\n"
         << "  \"sample_sha256\": \"" << sampleDigest << "\",\n"
         << "  \"normalization_source\": " << normalizationAudit << ",\n"
         << "  \"reaction\": \"gamma p -> p KShort KShort\",\n"
         << "  \"particle_order\": [\"beam\", \"recoil\", \"ks1\", "
            "\"ks2\"],\n"
         << "  \"polarization_fraction\": 0.0,\n"
         << "  \"normalization\": {\"generated_before_evaluation\": true, "
            "\"generated_events\": "
         << normalization->numGenEvents()
         << ", \"accepted_weight_sum\": " << normalization->numAccEvents()
         << ", \"s0_integral\": [" << norm.real() << ',' << norm.imag()
         << "]},\n"
         << "  \"event_count\": " << interface.numEvents() << ",\n"
         << "  \"weight_sum\": " << interface.sumWeights() << ",\n"
         << "  \"events\": [\n";

  reader->resetSource();
  for (int index = 0; index < interface.numEvents(); ++index) {
    std::unique_ptr<Kinematics> event(reader->getEvent());
    if (!event) {
      std::cerr << "data reader ended before the advertised event count\n";
      return 2;
    }
    const auto &particles = event->particleList();
    if (particles.size() != 4) {
      std::cerr << "event " << event->eventID()
                << " does not contain beam plus three final particles\n";
      return 2;
    }
    const TLorentzVector finalState =
        particles[1] + particles[2] + particles[3];
    const TLorentzVector pair = particles[2] + particles[3];
    output << "    {\"event_id\": " << event->eventID()
           << ", \"weight\": " << event->weight() << ", \"particles\": [";
    for (std::size_t particleIndex = 0; particleIndex < particles.size();
         ++particleIndex) {
      writeFourMomentum(output, particles[particleIndex]);
      if (particleIndex + 1 != particles.size())
        output << ',';
    }
    output << "], \"beam\": ";
    writeFourMomentum(output, particles[0]);
    output << ", \"final_state_sum\": ";
    writeFourMomentum(output, finalState);
    output << ", \"ksks_mass_gev\": " << pair.M()
           << ", \"intensity\": " << interface.intensity(index) << '}';
    output << (index + 1 == interface.numEvents() ? "\n" : ",\n");
  }
  output << "  ]\n}\n";
  std::cout << "wrote " << interface.numEvents()
            << " deterministic event intensities to " << outputPath << '\n';
  return 0;
}

void usage() {
  std::cerr << "usage:\n"
            << "  amptoolsReference prepare-norm CONFIG NORMINT AUDIT_JSON\n"
            << "  amptoolsReference evaluate CONFIG OUTPUT_JSON SAMPLE_SHA256 "
               "AUDIT_JSON\n";
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    usage();
    return 1;
  }
  const std::string command(argv[1]);
  if (command == "prepare-norm" && argc == 5)
    return prepareNormalization(argv[2], argv[3], argv[4]);
  if (command == "evaluate" && argc == 6) {
    const std::string sampleDigest(argv[4]);
    if (!isSha256(sampleDigest)) {
      std::cerr << "sample digest must be a 64-character SHA-256 value\n";
      return 1;
    }
    return evaluate(argv[2], argv[3], sampleDigest, argv[5]);
  }
  usage();
  return 1;
}
