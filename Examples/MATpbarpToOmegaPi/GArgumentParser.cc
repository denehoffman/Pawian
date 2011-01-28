/**
 * @file GArgumentParser.cpp
 */

/* Copyright (C) Dr. Ruediger Berlich and Karlsruhe Institute of Technology
 * (University of the State of Baden-Wuerttemberg and National Laboratory
 * of the Helmholtz Association)
 *
 * Contact: info [at] gemfony (dot) com
 *
 * This file is part of the Geneva library, Gemfony scientific's optimization
 * library.
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


#include "Examples/MATpbarpToOmegaPi/GArgumentParser.hh"
#include "ErrLogger/ErrLogger.hh"

namespace Gem
{
  namespace Geneva
  {
    /************************************************************************************************/
    /**
     * A function that parses the command line for all required parameters
     */
    bool parseCommandLine(
    		int argc, 
		char **argv,
		std::string& configFile
		  , boost::uint16_t& parallelizationMode
		  , bool& serverMode
		  , std::string& ip
		  , unsigned short& port
		  , Gem::Common::serializationMode& serMode
    ) {
      try{
		// Check the command line options. Uses the Boost program options library.
		po::options_description desc("Usage: evaluator [options]");
		desc.add_options()
		  ("help,h", "emit help message")
		  ("configFile,c", po::value<std::string>(&configFile)->default_value(DEFAULTCONFIGFILE),
		   "The name of the configuration file holding further configuration options")
		  ("parallelizationMode,p", po::value<boost::uint16_t>(&parallelizationMode)->default_value(DEFAULTPARALLELIZATIONMODE),
		   "Whether or not to run this optimization in serial mode (0), multi-threaded (1) or networked (2) mode")
		  ("serverMode,s","Whether to run networked execution in server or client mode. The option only gets evaluated if \"--parallelizationMode=2\"")
		  ("ip",po::value<std::string>(&ip)->default_value(DEFAULTIP), "The ip of the server")
		  ("port",po::value<unsigned short>(&port)->default_value(DEFAULTPORT), "The port of the server")
		  ("serMode", po::value<Gem::Common::serializationMode>(&serMode)->default_value(DEFAULTSERMODE),
		   "Specifies whether serialization shall be done in TEXTMODE (0), XMLMODE (1) or BINARYMODE (2)")
		  ;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		// Emit a help message, if necessary
		if (vm.count("help")) {
		  ErrMsg << desc;
		  return false;
		}

		serverMode=false;
		if (vm.count("parallelizationMode")) {
			if(parallelizationMode > 2) {
			  ErrMsg << "Error: the \"-p\" or \"--parallelizationMode\" option may only assume the\n"
				 << "values 0 (serial), 1 (multi-threaded) or 2 (networked). Leaving ..." << endmsg;
			  return false;
			}

			if(parallelizationMode == 2) if(vm.count("serverMode")) serverMode = true;
		}

		if(parallelizationMode != DEFAULTPARALLELIZATIONMODE ||  ip != DEFAULTIP  ||  port != DEFAULTPORT){
			std::string parModeString;
			switch(parallelizationMode) {
			case 0:
				parModeString = "serial";
				break;
			case 1:
				parModeString = "multi-threaded";
				break;
			case 2:
				parModeString = "networked";
				break;
			};

			Info << "\nRunning with the following command line options:\n"
			     << "configFile = " << configFile << "\n"
			     << "parallelizationMode = " << parModeString << "\n"
			     << "serverMode = " << (serverMode?"true":"false") << "\n"
			     << "ip = " << ip << "\n"
			     << "port = " << port << "\n"
			     << "serMode = " << serMode << "\n"
			     << endmsg;
		}
      }
      catch(...){
    	  ErrMsg << "Error parsing the command line" << endmsg;
    	  return false;
      }

      return true;
    }

    /************************************************************************************************/
    /**
     * A function that parses a config file for further parameters
     */
    bool parseConfigFile(
			 const std::string& configFile,
			 boost::uint16_t& nProducerThreads
			 , boost::uint16_t& nEvaluationThreads
			 , std::size_t& populationSize
			 , std::size_t& nParents
			 , boost::uint32_t& maxIterations
			 , long& maxMinutes
			 , boost::uint32_t& reportIteration
			 , recoScheme& rScheme
			 , sortingMode& smode
			 , std::size_t& arraySize
			 , boost::uint32_t& processingCycles
			 , bool& returnRegardless
			 , boost::uint32_t& waitFactor
			 , unsigned& jMax
			 , unsigned& pbarMom
			 , int& errLogMode

	) {
      boost::uint16_t recombinationScheme=0;
      bool verbose;

      // Check the name of the configuation file
      if(configFile.empty() || configFile == "empty" || configFile == "unknown") {
	ErrMsg << "Error: Invalid configuration file name given: \"" << configFile << "\"" << endmsg;
	return false;
      }

      try{
	// Check the configuration file line options. Uses the Boost program options library.
	po::options_description config("Allowed options");
	config.add_options()
	  ("nProducerThreads",po::value<boost::uint16_t>(&nProducerThreads)->default_value(DEFAULTNPRODUCERTHREADS),
	   "The amount of random number producer threads")
	  ("nEvaluationThreads",po::value<boost::uint16_t>(&nEvaluationThreads)->default_value(DEFAULTNEVALUATIONTHREADS),
	   "The amount of threads processing individuals simultaneously")
	  ("populationSize",po::value<std::size_t>(&populationSize)->default_value(DEFAULTPOPULATIONSIZE),
	   "The size of the super-population")
	  ("nParents",po::value<std::size_t>(&nParents)->default_value(DEFAULTNPARENTS),
	   "The number of parents in the population") // Needs to be treated separately
	  ("maxIterations", po::value<boost::uint32_t>(&maxIterations)->default_value(DEFAULTMAXITERATIONS),
	   "Maximum number of iterations in the population")
	  ("maxMinutes", po::value<long>(&maxMinutes)->default_value(DEFAULTMAXMINUTES),
	   "The maximum number of minutes the optimization of the population should run")
	  ("reportIteration",po::value<boost::uint32_t>(&reportIteration)->default_value(DEFAULTREPORTITERATION),
	   "The number of iterations after which information should be emitted in the super-population")
	  ("rScheme",po::value<boost::uint16_t>(&recombinationScheme)->default_value(DEFAULTRSCHEME),
	   "The recombination scheme for the super-population")
	  ("sortingScheme,o", po::value<sortingMode>(&smode)->default_value(DEFAULTSORTINGSCHEME),
	   "Determines whether sorting is done in MUCOMMANU (0), MUPLUSNU (1)  or MUNU1PRETAIN (2) mode")
	  ("arraySize", po::value<std::size_t>(&arraySize)->default_value(DEFAULTARRAYSIZE),
	   "The size of the buffer with random arrays in the random factory")
	  ("verbose",po::value<bool>(&verbose)->default_value(DEFAULTVERBOSE),
	   "Whether additional information should be emitted")
	  ("processingCycles", po::value<boost::uint32_t>(&processingCycles)->default_value(DEFAULTPROCESSINGCYCLES),
	   "The maximum number of cycles a client should perform mutations before it returns without success")
	  ("returnRegardless", po::value<bool>(&returnRegardless)->default_value(DEFAULTRETURNREGARDLESS),
	   "Specifies whether results should be returned even if they are not better than before")
	  ("waitFactor", po::value<boost::uint32_t>(&waitFactor)->default_value(DEFAULTGBTCWAITFACTOR),
	   "Influences the maximum waiting time of the GBrokerEA after the arrival of the first evaluated individuum")
	  ("jMax", po::value<unsigned>(&jMax)->default_value(DEFAULTJMAX),"choose jmax")
	  ("pbarMom", po::value<unsigned>(&pbarMom)->default_value(DEFAULTPBARMOM),"choose pbar momentum")
	  ("errLogMode", po::value<int>(&errLogMode)->default_value(DEFAULTERRLOGMODE),
	   "choose mode for Error logger")
	  ;
	
	po::variables_map vm;
	std::ifstream ifs(configFile.c_str());
	if(!ifs.good()) {
	  ErrMsg << "Error accessing configuration file " << configFile << endmsg;
	  return false;
	}

        store(po::parse_config_file(ifs, config), vm);
        notify(vm);

	// Emit a help message, if necessary
	if (vm.count("help")) {
	  Info << config << endmsg;
	  return false;
	}
	
	// Check the number of parents in the super-population
	if(2*nParents > populationSize){
	  ErrMsg << "Error: Invalid number of parents inpopulation\n"
		 << "nParents       = " << nParents << "\n"
		 << "populationSize = " << populationSize << endmsg;

	  return false;
	}

	// Workaround for assigment problem with rScheme
	if(recombinationScheme==(boost::uint16_t)VALUERECOMBINE)
	  rScheme=VALUERECOMBINE;
	else if(recombinationScheme==(boost::uint16_t)RANDOMRECOMBINE)
	  rScheme=RANDOMRECOMBINE;
	else if(recombinationScheme==(boost::uint16_t)DEFAULTRECOMBINE)
	  rScheme=DEFAULTRECOMBINE;
	else {
	  ErrMsg << "Error: Invalid recombination scheme in population: " << recombinationScheme << endmsg;
	  return false;
	}

	if(verbose){
	  Info << "\nRunning with the following options from " << configFile << ":\n"
	       << "nProducerThreads = " << (boost::uint16_t)nProducerThreads << "\n" // boost::uint8_t not printable on gcc ???
	       << "populationSize = " << populationSize << "\n"
	       << "nParents = " << nParents << "\n"
	       << "maxIterations = " << maxIterations << "\n"
	       << "maxMinutes = " << maxMinutes << "\n"
	       << "reportIteration = " << reportIteration << "\n"
	       << "rScheme = " << (boost::uint16_t)rScheme << "\n"
	       << "sortingScheme = " << smode << "\n"
	       << "arraySize = " << arraySize << "\n"
	       << "processingCycles = " << processingCycles << "\n"
	       << "returnRegardless = " << (returnRegardless?"true":"false") << "\n"
	       << "jMax = " << jMax << "\n"
	       << "pbarMom = " << pbarMom << "\n"
	       << "errLogMode = " << errLogMode << "\n"
	       << endmsg;
	}
      }
      catch(...){
	ErrMsg << "Error parsing the configuration file " << configFile << endmsg;
	return false;
      }

      return true;
    }

    /************************************************************************************************/

  } /* namespace GenEvA */
} /* namespace Gem */
