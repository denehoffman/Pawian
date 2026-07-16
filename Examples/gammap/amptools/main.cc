#include <iostream>

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/Kinematics.h"

#include "ROOTDataReader.h"
#include "Zlm.h"

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: printIntensity fit.cfg\n";
    return 1;
  }

  AmpToolsInterface::registerDataReader(ROOTDataReader());
  AmpToolsInterface::registerAmplitude(Zlm());

  ConfigFileParser parser(argv[1]);
  ConfigurationInfo *cfgInfo = parser.getConfigurationInfo();

  AmpToolsInterface ati(cfgInfo);

  unsigned int nEvents = ati.numEvents();

  for (unsigned int i = 0; i < nEvents; ++i) {
    Kinematics *kin = ati.kinematics(i);

    double intensity = ati.intensity(i);

    std::cout << i << " " << intensity << "\n";
  }

  return 0;
}
