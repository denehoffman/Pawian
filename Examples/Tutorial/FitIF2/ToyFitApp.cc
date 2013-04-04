//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>

// Boost header files go here
#include <boost/lexical_cast.hpp>

//#include "ErrLogger/ErrLogger.hh"

// Minimizer Interface header files go here
#include "Examples/Tutorial/FitIF2/MinimizerInterface/MIGeneric.hh"
#include "Examples/Tutorial/FitIF2/MinimizerInterface/MIGeneva.hh"
#include "Examples/Tutorial/FitIF2/MinimizerInterface/GArgumentParser.hh"

// Geneva header files go here
#include <hap/GRandomT.hpp>
#include <common/GCommonEnums.hpp>
#include <common/GExceptions.hpp>
#include <geneva/GConstrainedDoubleObject.hpp>
#include <geneva/GConstrainedDoubleObjectCollection.hpp>
#include <geneva/GDoubleGaussAdaptor.hpp>
#include <geneva/GObjectExpectationChecksT.hpp>
#include <geneva/GParameterObjectCollection.hpp>
#include <geneva/GParameterSet.hpp>
#include <courtier/GAsioHelperFunctions.hpp>
#include <courtier/GAsioTCPClientT.hpp>
#include <courtier/GAsioTCPConsumerT.hpp>
#include <geneva/GParameterSet.hpp>
#include <geneva/GBrokerEA.hpp>
#include <geneva/GEvolutionaryAlgorithm.hpp>
#include <geneva/GIndividual.hpp>
#include <geneva/GMultiThreadedEA.hpp>
#include <common/GSerializationHelperFunctionsT.hpp>
#include <geneva/GOptimizationEnums.hpp>

// The toy-data to fit to
#include "Examples/Tutorial/FitIF2/PolyFit.hh"

using namespace Gem::Geneva;
using namespace Gem::Courtier;
using namespace Gem::Hap;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  //std::string whichMinimizer("Geneva");
  double p0=-10., p1=10., p2=1., p3=-0.01, sigma_smear=3;

  // Generate data distribution
  boost::shared_ptr<MIGeneric> myFit(new PolyFit(p0, p1, p2, p3, sigma_smear));

  // Initiate parameters
  double val[4], min[4], max[4], err[4];
  val[0] = -11; max[0] = 0; min[0] = -20; err[0] = 3;
  val[1] = 9.8; max[1] = 15; min[1] = 5; err[1] = 2;
  val[2] = 1.1; max[2] = 1.5; min[2] = 0.5; err[2] = 0.3;
  val[3] = -0.008; max[3] = 0.; min[3] = -0.02; err[3] = 0.005; 
  myFit->setStartPar(4, val,  min, max, err);

  //--------------------------Minimizer IF --------------------------------------------------------
  //boost::shared_ptr<MIGeneva> minGen(new MIGeneva(myFit));

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

  configFile="./GStartProject.cfg"; parallelizationMode=1; serverMode=false; ip="localhost"; port=10000; serMode=Gem::Common::SERIALIZATIONMODE_TEXT;

  bool parsedConfig = parseConfigFile(configFile,
		      nProducerThreads,
		      nEvaluationThreads,
		      populationSize,
		      nParents,
		      maxIterations,
		      maxMinutes,
		      reportIteration,
		      rScheme,
		      smode,
		      arraySize,
		      processingCycles,
		      returnRegardless,
		      waitFactor);
    if(!parsedConfig) cout << "TODO" << endl;

////--------------------------------------------------------

  // Random numbers are our most valuable good. Set the number of threads
  GRANDOMFACTORY->setNProducerThreads(nProducerThreads);
  GRANDOMFACTORY->setArraySize(arraySize);

  //***************************************************************************
  // If this is a client in networked mode, we can just start the listener and
  // return when it has finished
  if(parallelizationMode==2 && !serverMode) {
    boost::shared_ptr<GAsioTCPClientT<GIndividual> > p(new GAsioTCPClientT<GIndividual>(ip, boost::lexical_cast<std::string>(port)));

    p->setMaxStalls(0); // An infinite number of stalled data retrievals
    p->setMaxConnectionAttempts(100); // Up to 100 failed connection attempts

    // Prevent return of unsuccessful adaption attempts to the server
    p->returnResultIfUnsuccessful(returnRegardless);

    // Start the actual processing loop
    p->run();

    return 0;
  }
  //***************************************************************************

  std::cout << "Create Parents ..." << std::endl << std::endl;

 // Create the first set of parent individuals. Initialization of parameters is (should be) done randomly.
  std::vector<boost::shared_ptr<MIGeneva> > parentIndividuals;
  for(std::size_t p = 0 ; p<nParents; p++) {
    //TODO: vary start parameters
    boost::shared_ptr<MIGeneva> gdii_ptr(new MIGeneva(myFit));
    gdii_ptr->setProcessingCycles(processingCycles);

    parentIndividuals.push_back(gdii_ptr);
  }

  std::cout << "Parent created ..." << std::endl << std::endl;

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
      boost::shared_ptr<GAsioTCPConsumerT<GIndividual> > gatc(new GAsioTCPConsumerT<GIndividual>(port));
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
  pop_ptr->setDefaultPopulationSize(populationSize,nParents);
  pop_ptr->setMaxIteration(maxIterations);
  pop_ptr->setMaxTime(boost::posix_time::minutes(maxMinutes));
  pop_ptr->setReportIteration(reportIteration);
  pop_ptr->setRecombinationMethod(rScheme);
  pop_ptr->setSortingScheme(smode);
  
  // Do the actual optimization
  pop_ptr->optimize();

  //--------------------------------------------------------------------------------------------

  //boost::shared_ptr<MIGeneva> bestIndividual_ptr=pop_ptr->getBestIndividual<MIGeneva>();


////--------------------------------------------------------

  //plot result
  //const MIGeneric& bestResult= bestIndividual_ptr->getResult();

  //std::cout << "Minimizer Geneva \t final par :\t" << bestResult.evaluator() << std::endl;
  //for(int i=0; i<bestResult.getNumPars(); i++){
  //  std::cout << "final Par" << i <<": \t" << bestResult.getPar(i) << " +- " << bestResult.getErr(i) << std::endl;
  //}

    /*std::cout << "Minimizer " << Nmin << "\t final par :\t" << genResult << std::endl;
    std::cout << "final a:\t" << val[0] << " +- " << err[0] << std::endl;
    std::cout << "final b:\t" << val[1] << " +- " << err[1] << std::endl; 
    std::cout << "final c:\t" << val[2] << " +- " << err[2] << std::endl;
    std::cout << "final d:\t" << val[3] << " +- " << err[3] << std::endl; 
    std::cout << "Done ..." << std::endl << std::endl;*/
  //}

  // Plot results
  //myFit->drawGraph(val[0],val[1],val[2],val[3]);
  return 0;
}
