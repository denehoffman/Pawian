#include <string>
#include <vector>

#include "GRandom.hpp"
#include "GEvolutionaryAlgorithm.hpp"
#include "GMultiThreadedEA.hpp"
#include "GBrokerEA.hpp"
#include "GIndividualBroker.hpp"
#include "GAsioTCPConsumer.hpp"
#include "GAsioTCPClient.hpp"
#include "GAsioHelperFunctions.hpp"

#include "Examples/MATpbarpToOmegaPi/OmegaPiData.hh"

// Declares a function to parse the command line
#include "Examples/MATpbarpToOmegaPi/GArgumentParser.hh"

// Information retrieval and printing
#include "Examples/MATpbarpToOmegaPi/GInfoFunction.hh"

#include "Examples/MATpbarpToOmegaPi/GOmegaPiIndividual.hh"
#include "Examples/MATpbarpToOmegaPi/PwaIntGeneva.hh"

using namespace std;

//TODO: 
//general fcn class

PwaIntGeneva::PwaIntGeneva(){
 // for(int i=in_paramNames.size()-1; i>-1; i--)
  //  paramNames.push_back(in_paramNames.at(i));
}

PwaIntGeneva::~PwaIntGeneva()
{
}

OmegaPiData::fitParamVal PwaIntGeneva::exec(std::vector<boost::shared_ptr<GOmegaPiIndividual> > parentIndividuals){
  
  //***************************************************************************
  // Create an instance of our optimization monitor, telling it to output information in given intervals
  std::ofstream resultSummary("./result.C");
  boost::shared_ptr<optimizationMonitor> om(new optimizationMonitor(nParents, resultSummary));

  //***************************************************************************
  // We can now start creating populations. We refer to them through the base class

  // This smart pointer will hold the different population types
  boost::shared_ptr<GEvolutionaryAlgorithm> pop_ptr;

  // Create the actual populations
  switch (parallelizationMode) {
    //-----------------------------------------------------------------------------------------------------
  case 0: // Serial execution
    // Create an empty population
    pop_ptr = boost::shared_ptr<GEvolutionaryAlgorithm>(new GEvolutionaryAlgorithm());
    break;

    //-----------------------------------------------------------------------------------------------------
  case 1: // Multi-threaded execution
    {
      // Create the multi-threaded population
      boost::shared_ptr<GMultiThreadedEA> popPar_ptr(new GMultiThreadedEA());

      // Population-specific settings
      popPar_ptr->setNThreads(nEvaluationThreads);

      // Assignment to the base pointer
      pop_ptr = popPar_ptr;
    }
    break;

    //-----------------------------------------------------------------------------------------------------
  case 2: // Networked execution (server-side)
    {
      // Create a network consumer and enrol it with the broker
      boost::shared_ptr<GAsioTCPConsumer> gatc(new GAsioTCPConsumer(port));
      gatc->setSerializationMode(serMode);
      GINDIVIDUALBROKER->enrol(gatc);

      // Create the actual broker population
      boost::shared_ptr<GBrokerEA> popBroker_ptr(new GBrokerEA());
      popBroker_ptr->setWaitFactor(waitFactor);

      // Assignment to the base pointer
      pop_ptr = popBroker_ptr;
    }
    break;
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Now we have suitable populations and can fill them with data

  // Add individuals to the population
  for(std::size_t p = 0 ; p<nParents; p++) {
    pop_ptr->push_back(parentIndividuals[p]);
  }
 
  // Specify some general population settings
  pop_ptr->setPopulationSize(populationSize,nParents);
  pop_ptr->setMaxIteration(maxIterations);
  pop_ptr->setMaxTime(boost::posix_time::minutes(maxMinutes));
  pop_ptr->setReportIteration(reportIteration);
  pop_ptr->setRecombinationMethod(rScheme);
  pop_ptr->setSortingScheme(smode);
  pop_ptr->registerInfoFunction(boost::bind(&optimizationMonitor::informationFunction, om, _1, _2));
  
  // Do the actual optimization
  pop_ptr->optimize();

  //--------------------------------------------------------------------------------------------

  // Make sure we close the result file
  resultSummary.close();

  boost::shared_ptr<GOmegaPiIndividual> bestIndividual_ptr=pop_ptr->getBestIndividual<GOmegaPiIndividual>();
  OmegaPiData::fitParamVal finalFitParm;
  assert(bestIndividual_ptr->getFitParams(finalFitParm));
  parentIndividuals.empty();
  parentIndividuals.push_back(bestIndividual_ptr);
 
  return finalFitParm;
}
