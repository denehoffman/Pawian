#ifndef _GPsi2STo2K2PiGamFit_H
#define _GPsi2STo2K2PiGamFit_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

// Boost header files go here
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>

#include "TROOT.h"
#include "Minuit2/MnUserParameters.h"

// GenEvA headers go here
#include "GEnums.hpp"
#include "GSerializationHelperFunctions.hpp"
#include "GEvolutionaryAlgorithm.hpp"

#include "Examples/Psi2STo2K2PiGam/GPsi2STo2K2PiGamIndividual.hh"
#include "Examples/Psi2STo2K2PiGam/GPsi2STo2K2PiGamInfoFunction.hh"

using namespace ROOT::Minuit2;
using namespace Gem::GenEvA;

class AbsPsi2STo2K2PiGamLh;

class GPsi2STo2K2PiGamFit{

public:

  // create/copy/destroy:

  ///Constructor 
  GPsi2STo2K2PiGamFit(boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLh,  MnUserParameters& upar);

  /** Destructor */
  virtual ~GPsi2STo2K2PiGamFit();

  void doFit();

  boost::shared_ptr<GPsi2STo2K2PiGamIndividual> bestIndividual() const {return _bestIndividual_ptr;};

protected:


private:

  boost::shared_ptr<AbsPsi2STo2K2PiGamLh> _lhPtr;
  MnUserParameters _upar;

  boost::uint16_t _parallelizationMode;
  bool _serverMode;
  std::string _ip;
  unsigned short _port;
  boost::uint16_t _nProducerThreads;
  boost::uint16_t _nEvaluationThreads;
  std::size_t _populationSize;
  std::size_t _nParents;
  boost::uint32_t _maxIterations;
  long _maxMinutes;
  boost::uint32_t _reportIteration;
  recoScheme _rScheme;
  std::size_t _arraySize;
  sortingMode _smode;
  boost::uint32_t _processingCycles;
  bool _returnRegardless;
  boost::uint32_t _waitFactor;
  std::vector<boost::shared_ptr<GPsi2STo2K2PiGamIndividual> > _parentIndividuals;
  boost::shared_ptr<optimizationMonitor> _optiMonitor;
  boost::shared_ptr<GEvolutionaryAlgorithm> _pop_ptr;
  std::ofstream _resultSummary;
  boost::shared_ptr<GPsi2STo2K2PiGamIndividual> _bestIndividual_ptr;
};

#endif
