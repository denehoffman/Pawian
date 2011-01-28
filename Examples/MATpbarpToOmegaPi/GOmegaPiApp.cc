/**
 * @file GStartProject.cpp
 */

/*
 * Copyright (C) Authors of the Geneva library collection and Karlsruhe
 * Institute of Technology (University of the State of Baden-Wuerttemberg
 * and National Laboratory of the Helmholtz Association).
 *
 * See the AUTHORS file in the top-level directory for a list of authors.
 *
 * Contact: info [at] gemfony (dot) com
 *
 * This file is part of the Geneva library collection
 *
 * Geneva is free software: you can redistribute it and/or modify
 * it under the terms of version 3 of the GNU Affero General Public License
 * as published by the Free Software Foundation.
 *
 * Geneva is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with the Geneva library.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For further information on Gemfony scientific and Geneva, visit
 * http://www.gemfony.com .
 */


// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>

//Root header files go here
#include <TStopwatch.h>

// Boost header files go here
#include <boost/lexical_cast.hpp>

// Geneva header files go here
#include <courtier/GAsioHelperFunctions.hpp>
#include <courtier/GAsioTCPClientT.hpp>
#include <courtier/GAsioTCPConsumerT.hpp>
#include <geneva/GBrokerEA.hpp>
#include <geneva/GEvolutionaryAlgorithm.hpp>
#include <geneva/GIndividual.hpp>
#include <geneva/GMultiThreadedEA.hpp>

#include "PwaUtils/pbarpStates.hh"
#include "Examples/MATpbarpToOmegaPi/pbarpToOmegaPi0States.hh"
// The individual that should be optimized
#include "Examples/MATpbarpToOmegaPi/GOmegaPiIndividual.hh"

// Declares a function to parse the command line
#include "Examples/MATpbarpToOmegaPi/GArgumentParser.hh"

// Information retrieval and printing
#include "Examples/MATpbarpToOmegaPi/GInfoFunction.hh"

#include "Examples/MATpbarpToOmegaPi/OmegaPiEventList.hh"
//#include "Examples/EtacToapi0Fit/EtacToapi0Lh.hh"
#include "Examples/MATpbarpToOmegaPi/OmegaPiHist.hh"

#include "Examples/MATpbarpToOmegaPi/OmegaPiData.hh"

#include "Setup/PwaEnv.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "Event/CBElsaReader.hh"
#include "Particle/PdtParser.hh"

#include "ErrLogger/ErrLogger.hh"

// The individual that should be optimized
//#include "GStartIndividual.hpp"

// Declares a function to parse the command line
//#include "GArgumentParser.hpp"

// Information retrieval and printing
//#include "GInfoFunction.hpp"

using namespace Gem::Geneva;
using namespace Gem::Courtier;
using namespace Gem::Hap;
//using namespace Gem::Util;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  std::string configFile;		  
  boost::uint16_t parallelizationMode;
  bool serverMode;
  std::string ip;
  unsigned short port;
  boost::uint16_t nProducerThreads;
  boost::uint16_t nEvaluationThreads;
  std::size_t populationSize;
  std::size_t nParents;
  boost::uint32_t maxIterations;
  long maxMinutes;
  boost::uint32_t reportIteration;
  recoScheme rScheme;
  std::size_t arraySize;
  sortingMode smode;
  boost::uint32_t processingCycles;
  bool returnRegardless;
  boost::uint32_t waitFactor;
  unsigned jMax;
  unsigned pbarMom;
  int errLogMode;

  Gem::Common::serializationMode serMode;

  TStopwatch timer;
  timer.Start();

  if(!parseCommandLine(argc, argv,
		       configFile,			  
		       parallelizationMode,
		       serverMode,
		       ip,
		       port,
		       serMode)
     ||
     !parseConfigFile(configFile,
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
		      waitFactor,
		      jMax,
		      pbarMom,
		      errLogMode))
    { exit(1); }

  // Random numbers are our most valuable good. Set the number of threads
  GRANDOMFACTORY->setNProducerThreads(nProducerThreads);
  GRANDOMFACTORY->setArraySize(arraySize);

  switch(errLogMode) {
  case -1:
    ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
    break;
  case 0:
    ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
    break;
  case 1:
    ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
    break;
  case 2:
    ErrLogger::instance()->setLevel(log4cpp::Priority::WARN);
    break;
  case 3:
    ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR);
    break;
  case 4:
    ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR);
    break;
  case 5:
    ErrLogger::instance()->setLevel(log4cpp::Priority::ALERT);
    break;
  default: 
    ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  }

  Info << "Maximum spin content: " << jMax << endmsg;
  Info << "pbar momentum: " << pbarMom   << endmsg;

  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
  
  std::string piomegaDatFile;
  std::string piomegaMcFile; 
  if(pbarMom==600){
    piomegaDatFile=theSourcePath+"/Examples/MATpbarpToOmegaPi/data/510_0600_sample.dat";
    piomegaMcFile=theSourcePath+"/Examples/MATpbarpToOmegaPi/data/mc510_0600_sample.dat";
  }
  else if (pbarMom==1940){
    piomegaDatFile=theSourcePath+"/Examples/MATpbarpToOmegaPi/data/510_1940.dat";
    piomegaMcFile=theSourcePath+"/Examples/MATpbarpToOmegaPi/data/mc510_1940.dat";
  }
  else{
    Alert <<"data for pbarMom= " << pbarMom << "not available; use 600 or 1940!!!" << endmsg;
    exit(1);
  }

  ParticleTable pTable;
  PdtParser parser;
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");
  if (!parser.parse(pdtFile, pTable)) {
    Alert << "Error: could not parse " << pdtFile << endmsg;
    exit(1);
  }

  std::vector<std::string> fileNames;

  fileNames.push_back(piomegaDatFile);
  CBElsaReader eventReader(fileNames, 3, 0); 
  EventList piOmegaEventsData;
  eventReader.fillAll(piOmegaEventsData);

  if (!piOmegaEventsData.findParticleTypes(pTable))
    Warning << "could not find all particles" << endmsg;

  Info << "\nFile has " << piOmegaEventsData.size() << " events. Each event has "
                  <<  piOmegaEventsData.nextEvent()->size() << " final state particles.\n" << endmsg;
  piOmegaEventsData.rewind();


  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = piOmegaEventsData.nextEvent()) != 0 && evtCount < 20) {
    Info << "\n" 
                    << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
                    << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
                    << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
                    << endmsg;
    ++evtCount;
  }
  piOmegaEventsData.rewind();

  std::vector<std::string> fileNamesMc;

  fileNamesMc.push_back(piomegaMcFile);
  CBElsaReader eventReaderMc(fileNamesMc, 3, 0); 
  EventList piOmegaEventsMc;
  eventReaderMc.fillAll(piOmegaEventsMc);
  piOmegaEventsMc.rewind();

  //boost::shared_ptr<const OmegaPiEventList> theOmegaPiEventPtr(new OmegaPiEventList(piOmegaEventsData, piOmegaEventsMc, jMax,  pbarMom));
  boost::shared_ptr<OmegaPiEventList> theOmegaPiEventPtr = OmegaPiEventList::getList(piOmegaEventsData, piOmegaEventsMc, jMax,  pbarMom);
  //theOmegaPiEventPtr->initList(piOmegaEventsData, piOmegaEventsMc, jMax,  pbarMom);

  boost::shared_ptr<pbarpStates> pbarpStatesPtr(new pbarpStates(jMax));
  //boost::shared_ptr<pbarpToOmegaPi0States> pbarpToOmegaPi0StatesPtr(new pbarpToOmegaPi0States(pbarpStatesPtr));
  boost::shared_ptr<pbarpToOmegaPi0States> pbarpToOmegaPi0StatesPtr = pbarpToOmegaPi0States::getStates(pbarpStatesPtr);
  
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

  // Create the first set of parent individuals. Initialization of parameters is done randomly.
  std::vector<boost::shared_ptr<GOmegaPiIndividual> > parentIndividuals;
  for(std::size_t p = 0 ; p<nParents; p++) {
    boost::shared_ptr<GOmegaPiIndividual> gdii_ptr(new GOmegaPiIndividual(pbarpToOmegaPi0StatesPtr));
    gdii_ptr->setProcessingCycles(processingCycles);

    parentIndividuals.push_back(gdii_ptr);
  }

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
  pop_ptr->registerInfoFunction(boost::bind(&optimizationMonitor::informationFunction, om, _1, _2));
  
  // Do the actual optimization
  pop_ptr->optimize();

  //--------------------------------------------------------------------------------------------

  // Make sure we close the result file
  resultSummary.close();

  boost::shared_ptr<GOmegaPiIndividual> bestIndividual_ptr=pop_ptr->getBestIndividual<GOmegaPiIndividual>();
  OmegaPiData::fitParamVal finalFitParm;
  assert(bestIndividual_ptr->getFitParams(finalFitParm));
  boost::shared_ptr<OmegaPiLh> finalOmegaPiLh=bestIndividual_ptr->getOmegaPiLhPtr();
  OmegaPiHist theHistogrammer(finalOmegaPiLh,finalFitParm);

  Info << "Final LH:\n\n";

  bestIndividual_ptr->printFitParams(finalFitParm);

  Info << "Done ...\n";

  timer.Stop();
  Info << "Real time: " << (double)timer.RealTime() << " s, CPU time: " << (double)timer.CpuTime() << " s\a\a\n";

  return 0;
}
