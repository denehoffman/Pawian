#ifndef _MIGENEVA_H
#define _MIGENEVA_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include "Examples/Tutorial/FitIF/MinimizerInterface/MIData.hh"
#include "Examples/Tutorial/FitIF/MinimizerInterface/MIBase.hh"
#include "Examples/Tutorial/FitIF/MinimizerInterface/GArgumentParser.hh"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// Geneva header files go here
#include <courtier/GAsioHelperFunctions.hpp>
#include <courtier/GAsioTCPClientT.hpp>
#include <courtier/GAsioTCPConsumerT.hpp>
#include <geneva/GBrokerEA.hpp>
#include <geneva/GEvolutionaryAlgorithm.hpp>
#include <geneva/GIndividual.hpp>
#include <geneva/GMultiThreadedEA.hpp>
#include <common/GCommonEnums.hpp>
#include <common/GSerializationHelperFunctionsT.hpp>
#include <geneva/GOptimizationEnums.hpp>

using namespace std;

class MIGeneva : public MIBase {

public:
  /// Default Constructor (0x0)
  MIGeneva(boost::shared_ptr<MIData> theData, string inConfigFile="./GStartProject.cfg", boost::uint16_t inparallelizationMode=1, bool inserverMode=false, string inip="localhost", unsigned short inport=10000, Gem::Common::serializationMode inserMode=Gem::Common::SERIALIZATIONMODE_TEXT);  
  virtual const double exec(unsigned int Npar, double* par,  double* min, double* max, double* err); 

  /** Destructor */
  virtual ~MIGeneva();

 protected:

 private:
   boost::shared_ptr<MIData> _myData;
 // vector<string> paramNames;

  std::string configFile;	
  boost::uint16_t parallelizationMode;
  bool serverMode;
  std::string ip;
  unsigned short port;
  Gem::Common::serializationMode serMode;
  boost::uint16_t nProducerThreads;
  boost::uint16_t nEvaluationThreads;
  std::size_t populationSize;
  std::size_t nParents;
  boost::uint32_t maxIterations;
  long maxMinutes;
  boost::uint32_t reportIteration;
  Gem::Geneva::recoScheme rScheme;
  std::size_t arraySize;
  Gem::Geneva::sortingMode smode;
  boost::uint32_t processingCycles;
  bool returnRegardless;
  boost::uint32_t waitFactor;

};

#endif /* _MIGeneva_H */
