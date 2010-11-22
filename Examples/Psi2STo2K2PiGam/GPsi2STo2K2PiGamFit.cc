#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/GPsi2STo2K2PiGamFit.hh"
#include "Examples/Psi2STo2K2PiGam/AbsPsi2STo2K2PiGamLh.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Minuit2/MnUserParameters.h"

#include "GMultiThreadedEA.hpp"
#include "GBrokerEA.hpp"
#include "GIndividualBroker.hpp"
#include "GAsioTCPConsumer.hpp"
#include "GAsioTCPClient.hpp"
#include "GAsioHelperFunctions.hpp"

GPsi2STo2K2PiGamFit::GPsi2STo2K2PiGamFit(boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLh, MnUserParameters& upar) :
  _lhPtr(theLh),
  _upar(upar),
  _parallelizationMode(1),
  _nProducerThreads(8),
  _nEvaluationThreads(16),
  _populationSize(50),
  _nParents(2),
  _maxIterations(100),
  _maxMinutes(0),
  _reportIteration(1),
  _rScheme(VALUERECOMBINE),
  _arraySize(1000),
  _smode(MUCOMMANU),
  _processingCycles(1),
  _returnRegardless(0),
  _waitFactor(2),
  _resultSummary("./result.C")
{
  GRANDOMFACTORY->setNProducerThreads(_nProducerThreads);
  GRANDOMFACTORY->setArraySize(_arraySize);

  for(std::size_t p = 0 ; p<_nParents; p++) {
    boost::shared_ptr<GPsi2STo2K2PiGamIndividual> gdii_ptr(new GPsi2STo2K2PiGamIndividual(theLh, upar));
    gdii_ptr->setProcessingCycles(_processingCycles);
    
    _parentIndividuals.push_back(gdii_ptr);
  }

  _optiMonitor=boost::shared_ptr<optimizationMonitor>(new optimizationMonitor(_nParents, _resultSummary));


  // Create the multi-threaded population
  boost::shared_ptr<GMultiThreadedEA> popPar_ptr(new GMultiThreadedEA());
  
  // Population-specific settings
  popPar_ptr->setNThreads(_nEvaluationThreads);
  
  // Assignment to the base pointer
  _pop_ptr = popPar_ptr;


  // Now we have suitable populations and can fill them with data

  // Add individuals to the population
  for(std::size_t p = 0 ; p<_nParents; p++) {
    _pop_ptr->push_back(_parentIndividuals[p]);
  }
 
  // Specify some general population settings
  _pop_ptr->setPopulationSize(_populationSize,_nParents);
  _pop_ptr->setMaxIteration(_maxIterations);
  _pop_ptr->setMaxTime(boost::posix_time::minutes(_maxMinutes));
  _pop_ptr->setReportIteration(_reportIteration);
  _pop_ptr->setRecombinationMethod(_rScheme);
  _pop_ptr->setSortingScheme(_smode);
  _pop_ptr->registerInfoFunction(boost::bind(&optimizationMonitor::informationFunction,_optiMonitor, _1, _2));
}


GPsi2STo2K2PiGamFit::~GPsi2STo2K2PiGamFit()
{
}


void GPsi2STo2K2PiGamFit::doFit(){ 
  _pop_ptr->optimize();

  // Make sure we close the result file
  _resultSummary.close();
  _bestIndividual_ptr=_pop_ptr->getBestIndividual<GPsi2STo2K2PiGamIndividual>();
}
