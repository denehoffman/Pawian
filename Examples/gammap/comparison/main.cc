#include <iostream>

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "MinuitInterface/MinuitMinimizationManager.h"

#include "ROOTDataReader.h"
#include "Zlm.h"

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: fitZlm fit.cfg\n";
    return 1;
  }

  AmpToolsInterface::registerDataReader(ROOTDataReader());
  AmpToolsInterface::registerAmplitude(Zlm());

  ConfigFileParser parser(argv[1]);
  ConfigurationInfo *cfgInfo = parser.getConfigurationInfo();

  AmpToolsInterface ati(cfgInfo);
  MinuitMinimizationManager *minuit = ati.minuitMinimizationManager();
  minuit->migradMinimization();
  ati.finalizeFit();
  return minuit->status() == 0 ? 0 : 2;
}
