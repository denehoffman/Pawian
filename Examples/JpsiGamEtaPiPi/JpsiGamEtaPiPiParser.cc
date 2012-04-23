// Bertram Kopf (RUB)

#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;
  
/************************************************************************************************/
/**
  * A function that parses the command line for all required parameters
  */
bool JpsiGamEtaPiPiParser::parseCommandLine(int argc, char **argv)
{
  try
  {
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
      ("configFile,c",po::value<std::string>(&_configFile)->default_value(_configFile),
	    "The name of the configuration file holding further configuration options")
      ;

    po::options_description common("Common Options");
    common.add_options()
      ("errLogMode,e", po::value<string>(&strErrLogMode)->default_value(strErrLogMode),"choose mode for Error logger.")
      ("datFile",po::value<string>(&_dataFile), "full path of data file")
      ("mcFile",po::value<string>(&_mcFile), "full path of Monte Carlo file")
      ("paramFile",po::value<string>(&_paramFile), "file with start parameters for fit or QA (full path)")
      ("startHypo",po::value<string>(&_startHypo), "choose the hyopthesis to start")
      ("enableHyp",po::value< vector<string> >(&_enabledHyps), "enable hypotheses")
      ("mode",po::value<string>(&_mode), "modes are: pwa, dumpDefaulParams, qaMode, plotmode")
      ("massIndependentFit", po::value<bool>(&_massIndependentFit), "enable/disable mass independence in fit")
      ("commonProdPhases",po::value<bool>(&_useCommonProductionPhases), "enable/disable common production phases")
      ;
    
    po::options_description config("Configuration file options");
    config.add_options()
      ("verbose",po::value<bool>(&verbose)->default_value(true), "Determines whether additional information should be emitted")
      ("mnParFix",po::value< vector<string> >(&_mnParFixs),  "minuit parameters can be fixed here")
      ("massRangeMin",po::value<double>(&_massMin), "min of eta pi pi mass range for mass indep. fit")
      ("massRangeMax",po::value<double>(&_massMax), "max of eta pi pi mass range for mass indep. fit")
      ;
    
    po::options_description cmdline_options;
    cmdline_options.add(desc).add(common);

    po::options_description config_file_options;
    config_file_options.add(config).add(common);
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
    po::notify(vm);

    // Check the name of the configuation file
    if(_configFile.empty() || _configFile == "empty" || _configFile == "unknown") 
    {
      std::cout << cmdline_options << endl;
      stringstream strError;
      strError << "Error: Invalid configuration file name given: \"" << _configFile << "\"";
      throw runtime_error(strError.str());
    }

    std::ifstream ifs(_configFile.c_str());
    if(!ifs.good()) 
    {
      stringstream strError;
      strError << "Error accessing configuratiocommonn file " << _configFile;
      std::cout << cmdline_options << endl;

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

	
    if(strErrLogMode == "debug") _errLogMode = debug;
    else if(strErrLogMode == "trace") _errLogMode = trace;
    else if(strErrLogMode == "routine") _errLogMode = routine;
    else if(strErrLogMode == "warning") _errLogMode = warning;
    else if(strErrLogMode == "error")   _errLogMode = error; 
    else if(strErrLogMode == "alert")   _errLogMode = alert; 
    else 
    {
      _errLogMode = debug;
      Warning << "ErrorLogger not (properly) set -> Use mode 'DEBUG' " ;  // << endmsg;
    }

    
    if(verbose){
      std::cout << "\nRunning with the following options using " << _configFile << ":\n\n"
		<< "Error log mode: " << _errLogMode <<"\n\n"
		<< "data file: " << _dataFile <<"\n\n"
		<< "mc file: " << _mcFile <<"\n\n"
		<< "file with start parameters for fit or qa: " << _paramFile << "\n\n"
		<< "startHypo: " << _startHypo << "\n\n"
		<< "mode: " << _mode << "\n\n"
	    << endl;


      std::vector<std::string>::const_iterator it;

      for (it=_enabledHyps.begin(); it!=_enabledHyps.end();++it){
    	  std::cout << "hypothesis\t" << (*it) << "\t enabled\n";
      }
      std::cout << std::endl;
    

      for (it=_mnParFixs.begin(); it!=_mnParFixs.end();++it){
    	  std::cout << "minuit parameter\t" << (*it) << "\t fixed\n";
      }
      std::cout << std::endl;
    }
  }  


  catch( std::exception & e )
  {
    cerr << "Error parsing the command line:" << endl;
    cerr << e.what() << std::endl;
    cerr << "You can use -h or --help to obtain the description of the program parameters." << endl;
    cerr << "This is the command line options\n" << endl;
    return false;
  }
  catch(...){
      std::cerr << "Error parsing the command line. Use -h or --help to see the description of the program paramters." << endl;
      return false;
  }

  return true;
}
