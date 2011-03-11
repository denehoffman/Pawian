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


#include "Examples/pbarpToOmegaPi/GArgumentParser.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace std;
  
/************************************************************************************************/
/**
  * A function that parses the command line for all required parameters
  */
bool ApplicationParameter::parseCommandLine(int argc, char **argv)
{
  try
  {
    int iAppMode=0;
    unsigned short DEFAULTPORT=port;
    std::string DEFAULTIP=ip;
    bool DEFAULTPARALLELIZATIONMODE=parallelizationMode;
    boost::uint16_t recombinationScheme=0;
    bool verbose;
    string strErrLogMode="debug";

    // Check the command line options. Uses the Boost program options library.
    string strAppName(argv[0]);
    size_t found = strAppName.rfind("/")+1;
    if (found != string::npos) strAppName=strAppName.substr(found);
    
    string strDesc="Usage: " + strAppName + " [options]";
    po::options_description desc(strDesc);
    desc.add_options()
      ("help,h", "emit help message")
      ("configFile,c",po::value<std::string>(&configFile)->default_value(configFile),
	    "The name of the configuration file holding further configuration options")
      ("parallelizationMode,p", po::value<boost::uint16_t>(&parallelizationMode)->default_value(parallelizationMode),
	    "Whether or not to run this optimization in serial mode (0), multi-threaded (1) or networked (2) mode")
      ("serverMode,s","Whether to run networked execution in server or client mode. The option only gets evaluated if \"--parallelizationMode=2\"")
      ("ip",po::value<std::string>(&ip)->default_value(ip), "The ip of the server")
      ("port",po::value<unsigned short>(&port)->default_value(port), "The port of the server")
      ("serMode", po::value<Gem::Common::serializationMode>(&serMode)->default_value(serMode),
	    "Specifies whether serialization shall be done in TEXTMODE (0), XMLMODE (1) or BINARYMODE (2)")
      ("spf", po::value<string>(&strPathStartParamFile),"Path to the start parameter file.")
      ("M1", "Spin M for calculation of spin density.")
      ("M2", "Spin M' for calculation of spin density.")
      ("allSpin","Calculate all possible spin density elements of the spin density matrix.")
      ("remB",po::value<int>(&nNumEventsRemove),"Remove number of events from the beginning of the event list.")
      ("remE",po::value<int>(&nNumEventsRemove),"Remove number of events from the end of the event list.")
      ("redE",po::value<int>(&nNumEventsRed),"Reduce number of events from the event list to specified number.")
      ;

    po::options_description common("Common Options");
    common.add_options()
      ("jmax", po::value<unsigned>(&jMax)->default_value(jMax),"choose jmax.")
      ("pbarmom", po::value<unsigned>(&pbarMom)->default_value(pbarMom),"choose pbar momentum.")
      ("errLogMode,e", po::value<string>(&strErrLogMode)->default_value(strErrLogMode),
      "choose mode for Error logger.")	  
      ("appMode,a", po::value<int>(&iAppMode)->default_value(enAppExecMode), 
       "Specifies which minimization algorithm should be used for fitting GenEvA (0) or Minuit (1) or Minuit after GenEvA (2) or spin density calculation (3) or QA mode (4).")
      ("name,n", po::value<string>(&strName),
        "Name that is attached to all otuput file names to be able to run multiple fits in parallel.")
      ;

    po::options_description config("Configuration file options");
    config.add_options()
      ("nProducerThreads",po::value<boost::uint16_t>(&nProducerThreads)->default_value(nProducerThreads),
      "The amount of random number producer threads")
      ("nEvaluationThreads",po::value<boost::uint16_t>(&nEvaluationThreads)->default_value(nEvaluationThreads),
      "The amount of threads processing individuals simultaneously")
      ("populationSize",po::value<std::size_t>(&populationSize)->default_value(populationSize),
      "The size of the super-population")
      ("nParents",po::value<std::size_t>(&nParents)->default_value(nParents),
      "The number of parents in the population") // Needs to be treated separately
      ("maxIterations", po::value<boost::uint32_t>(&maxIterations)->default_value(maxIterations),
      "Maximum number of iterations in the population")
      ("maxMinutes", po::value<long>(&maxMinutes)->default_value(maxMinutes),
      "The maximum number of minutes the optimization of the population should run")
      ("reportIteration",po::value<boost::uint32_t>(&reportIteration)->default_value(reportIteration),
      "The number of iterations after which information should be emitted in the super-population")
      ("rScheme",po::value<boost::uint16_t>(&recombinationScheme)->default_value(VALUERECOMBINE),
      "The recombination scheme for the super-population")
      ("sortingScheme,o", po::value<sortingMode>(&smode)->default_value(smode),
      "Determines whether sorting is done in MUCOMMANU (0), MUPLUSNU (1)  or MUNU1PRETAIN (2) mode")
      ("arraySize", po::value<std::size_t>(&arraySize)->default_value(arraySize),
      "The size of the buffer with random arrays in the random factory")
      ("verbose",po::value<bool>(&verbose)->default_value(true),
      "Whether additional information should be emitted")
      ("processingCycles", po::value<boost::uint32_t>(&processingCycles)->default_value(processingCycles),
      "The maximum number of cycles a client should perform mutations before it returns without success")
      ("returnRegardless", po::value<bool>(&returnRegardless)->default_value(returnRegardless),
      "Specifies whether results should be returned even if they are not better than before")
      ("waitFactor", po::value<boost::uint32_t>(&waitFactor)->default_value(waitFactor),
      "Influences the maximum waiting time of the GBrokerEA after the arrival of the first evaluated individuum")
      ("SourcePath", po::value<std::string>(&theSourcePath)->default_value(theSourcePath),
       "Specifies the path to root directory of the source")
      ;

    po::options_description cmdline_options;
    cmdline_options.add(desc).add(common);

    po::options_description config_file_options;
    config_file_options.add(config).add(common);
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
    po::notify(vm);

    // Check the name of the configuation file
    if(configFile.empty() || configFile == "empty" || configFile == "unknown") 
    {
      stringstream strError;
      strError << "Error: Invalid configuration file name given: \"" << configFile << "\"";
      throw runtime_error(strError.str());
    }

    std::ifstream ifs(configFile.c_str());
    if(!ifs.good()) 
    {
      stringstream strError;
      strError << "Error accessing configuratiocommonn file " << configFile;
      throw runtime_error(strError.str());
    }

    store(po::parse_config_file(ifs, config_file_options), vm);
    po::notify(vm);

    // Emit a help message, if necessary
    if (vm.count("help")) 
    {
      std::cout << cmdline_options << endl;
      exit(0);
    }    

    serverMode=false;
    if (vm.count("parallelizationMode")) {
	    if(parallelizationMode > 2) {
         stringstream strError;
         strError << "Error: the \"-p\" or \"--parallelizationMode\" option may only assume the\n"
                  << "values 0 (serial), 1 (multi-threaded) or 2 (networked). Leaving ...";
         throw runtime_error(strError.str());
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

	    std::cout << "\nRunning with the following command line options:\n"
		      << "configFile = " << configFile << "\n"
		      << "parallelizationMode = " << parModeString << "\n"
		      << "serverMode = " << (serverMode?"true":"false") << "\n"
		      << "ip = " << ip << "\n"
		      << "port = " << port << "\n"
		      << "serMode = " << serMode << "\n"
		      << endl;
    }
	
    enAppExecMode=enExecMode(iAppMode);
    if (vm.count("appMode"))
    {
      if (int(enAppExecMode) > 4 || int(enAppExecMode) < 0)
      {
        stringstream strError;
        strError << "The option number " << int(enAppExecMode) << " for minimization algorithm is out of range!";
        throw runtime_error(strError.str());
      }
    }
    
    if (vm.count("allSpin"))
    {
      bCalcAllSpindensity = true;
    }
    else bCalcAllSpindensity = false;
    
   if (vm.count("spf"))
   {
     if (strPathStartParamFile.empty()) 
     {
       stringstream strError;
       strError << "start parameter file must be specified when using --spf";
       throw runtime_error(strError.str());
     }
   }

   if (!vm.count("name") && (enAppExecMode != SpinDensity))
   {
    stringstream strError;
    strError << "Name of the process must be specified!";
    throw runtime_error(strError.str());
   }

  //Config File Parameter check
    
    // Check the number of parents in the super-population
    if(2*nParents > populationSize){
      stringstream strError;
      strError << "Error: Invalid number of parents inpopulation\n"
          << "nParents       = " << nParents << "\n"
          << "populationSize = " << populationSize;
      throw runtime_error(strError.str());
    }

    // Workaround for assigment problem with rScheme
    if(recombinationScheme==(boost::uint16_t)VALUERECOMBINE)
      rScheme=VALUERECOMBINE;
    else if(recombinationScheme==(boost::uint16_t)RANDOMRECOMBINE)
      rScheme=RANDOMRECOMBINE;
    else if(recombinationScheme==(boost::uint16_t)DEFAULTRECOMBINE)
      rScheme=DEFAULTRECOMBINE;
    else {
      stringstream strError;
      strError << "Error: Invalid recombination scheme in population: " << recombinationScheme;
      throw runtime_error(strError.str());
    }

    if(strErrLogMode == "debug") errLogMode = debug;
    else if(strErrLogMode == "trace") errLogMode = trace;
    else if(strErrLogMode == "routine") errLogMode = routine;
    else if(strErrLogMode == "warning") errLogMode = warning;
    else if(strErrLogMode == "error")   errLogMode = error; 
    else if(strErrLogMode == "alert")   errLogMode = alert; 
    else 
    {
      errLogMode = debug;
      Warning << "ErrorLogger not (properly) set -> Use mode 'DEBUG' " << endmsg;
    }

    if(vm.count("remB") && vm.count("remE")) 
    {
      stringstream strError;
      strError << "Options remB and remE can't be combined only one can be used.";
      throw runtime_error(strError.str());
    }
    else
    {
      if(vm.count("remE"))
      {
        bRemoveEventsFromEnd = true;
      }
      else if(vm.count("remB"))
      {
        bRemoveEventsFromEnd = false;
      }
      
      if (nNumEventsRemove<0)
      {
        stringstream strError;
        strError << "Number of events to remove must be positive!";
        throw runtime_error(strError.str());
      }
    }
    
    if (vm.count("redE"))
    {
      if(vm.count("remB") || vm.count("remE")) 
      {
        stringstream strError;
        strError << "Options redE can't be used in connection with options remE or remB.";
        throw runtime_error(strError.str());
      }
      
      if (nNumEventsRed < 0)
      {
        stringstream strError;
        strError << "Number to which number of events should be reduced must be greater then zero.";
        throw runtime_error(strError.str());
      }
      
      bRemoveEventsFromEnd = true;
    }
    else nNumEventsRed=0;

    if (vm.count("allSpin") || vm.count("M1") || vm.count("M2"))
    {
      if (!vm.count("spf"))
      {
        stringstream strError;
        strError << "Start parameter file must be specified for the calculation of spin density";
        throw runtime_error(strError.str());
      }
    }
    
    
    if(verbose){
      std::cout << "\nRunning with the following options using " << configFile << ":\n"
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
	    << "errLogMode = " << strErrLogMode << "\n"
	    << "SourcePath = " << theSourcePath << "\n"
	    << endl;
    }

  }

  catch( std::exception & e )
  {
    cerr << "Error parsing the command line:" << endl;
    cerr << e.what() << std::endl;
    cerr << "You can use -h or --help to obtain the description of the program parameters." << endl;
    return false;
  }
  catch(...){
      std::cerr << "Error parsing the command line. Use -h or --help to see the description of the program paramters." << endl;
      return false;
  }

  return true;
}
