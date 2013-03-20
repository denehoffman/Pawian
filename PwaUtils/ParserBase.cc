// Bertram Kopf (RUB)

#include "PwaUtils/ParserBase.hh"
#include "ErrLogger/ErrLogger.hh"
#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;


ParserBase::ParserBase(int argc,char **argv)
      : _configFile("./JpsiGamEtaPiPi.cfg")
      , _errLogMode(debug)
      , _dataFile("")
      , _mcFile("")
      , _paramFile("")
      , _serializationFile("")
      , _startHypo("base")
      , _mode("plotmode")
      , _outputFileNameSuffix("")
      , _verbose(true)
      , _noOfThreads(16)
      , _ratioMcToData(100000)
      , _cacheAmps(true)
      , _strErrLogMode("debug")
      , _desc(0)
      , _common(new po::options_description("Common Options"))
      , _config(new po::options_description("Configuration file options"))
      ,_useEvtWeight(false)
      ,_usePhaseSpaceHyp(false)
      ,_pdgTableFile("/Particle/pdtNew.table")
      ,_productionFormalism("Cano")
     {
       string strErrLogMode="debug";
    // Check the command line options. Uses the Boost program options library.
    string strAppName(argv[0]);
    size_t found = strAppName.rfind("/")+1;
    if (found != string::npos) strAppName=strAppName.substr(found);
    
    string strDesc="Usage: " + strAppName + " [options]";
    _desc= new po::options_description(strDesc);

    _desc->add_options()
      ("help,h", "emit help message")
      ("configFile,c",po::value<std::string>(&_configFile)->default_value(_configFile),
	    "The name of the configuration file holding further configuration options")
      ;

    _common->add_options()
      ("errLogMode,e", po::value<string>(&_strErrLogMode)->default_value(_strErrLogMode),"choose mode for Error logger.")
      ("datFile",po::value<string>(&_dataFile), "full path of data file")
      ("mcFile",po::value<string>(&_mcFile), "full path of Monte Carlo file")
      ("paramFile",po::value<string>(&_paramFile), "file with start parameters for fit or QA (full path)")
      ("serializationFile", po::value<string>(&_serializationFile), "serialized pwa i/o file")
      ("startHypo",po::value<string>(&_startHypo), "choose the hyopthesis to start")
      ("mode",po::value<string>(&_mode), "modes are: pwa, dumpDefaultParams, qaMode, plotmode, spinDensity")
      ("noOfThreads",po::value<int>(&_noOfThreads),  "number of threads for multi threaded mode")
      ("ratioMcToData",po::value<int>(&_ratioMcToData),  "number of MC events defined by ratio #MCs/#Data")
      ("cacheAmps",po::value<bool>(&_cacheAmps),  "cache amplitudes")
      ("useEventWeight",po::value<bool>(&_useEvtWeight), "enable/disable input for event weight")
      ("usePhaseSpaceHyp",po::value<bool>(&_usePhaseSpaceHyp), "use hypothesis for phase space") 
      ("name",po::value<string>(&_outputFileNameSuffix), "name that is attached to all otuput file names")
      ("pdgTableFile",po::value<string>(&_pdgTableFile), "path of the pdg-table file relative to the top dir")
      ;

    _config->add_options()
      ("verbose",po::value<bool>(&_verbose)->default_value(true), "Determines whether additional information should be emitted")
      ("enableHyp",po::value< vector<string> >(&_enabledHyps), "enable hypotheses")
      ("mnParFix",po::value< vector<string> >(&_mnParFixs),  "minuit parameters can be fixed here")
      ("finalStateParticle",po::value< vector<string> >(&_finalStateParticles),  "name of final state particles")
      ("decay",po::value< vector<string> >(&_decaySystem),  "decay: mother and pair of decay particles")
      ("addDynamics",po::value< vector<string> >(&_dynamics), "add dynamics/line shape for resonances")
      ("replaceParamSuffix",po::value< vector<string> >(&_replaceParSuffix),  "replace suffix for fit parameter name")
      ("replaceMassKey",po::value< vector<string> >(&_replaceMassKey),  "replace Key for the fit parameter of the mass")
      ("production",po::value< vector<string> >(&_productionSystem),  "pair of produced particles")
      ("productionFormalism",po::value< string >(&_productionFormalism),  "used formalism for the production")
      ("cloneParticle",po::value< vector<string> >(&_cloneParticle),  "particles to be cloned")
      ("histMass",po::value< vector<string> >(&_histMass),  "histograms inv mass for the selected final state paricles")
      ("histAngles",po::value< vector<string> >(&_histAngles),  "histograms decay angles")
      ;



  }  

/************************************************************************************************/
/************************************************************************************************/
/**
  * A function that parses the command line for all required parameters
  */

bool ParserBase::parseCommandLine(int argc, char **argv)
{
  try
  {
    
    po::options_description cmdline_options;
    cmdline_options.add(*_desc).add(*_common);

    po::options_description config_file_options;
    config_file_options.add(*_config).add(*_common);
    
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

        
    if(_strErrLogMode == "debug") _errLogMode = debug;
    else if(_strErrLogMode == "trace") _errLogMode = trace;
    else if(_strErrLogMode == "routine") _errLogMode = routine;
    else if(_strErrLogMode == "warning") _errLogMode = warning;
    else if(_strErrLogMode == "error")   _errLogMode = error; 
    else if(_strErrLogMode == "alert")   _errLogMode = alert; 
    else 
    {
      _errLogMode = debug;
      Warning << "ErrorLogger not (properly) set -> Use mode 'DEBUG' " ;  // << endmsg;
    }

    
    if(_verbose){
      std::cout << "\nRunning with the following options using " << _configFile << ":\n\n"
                << "Error log mode: " << _errLogMode <<"\n\n"
                << "data file: " << _dataFile <<"\n\n"
                << "mc file: " << _mcFile <<"\n\n"
                << "file with start parameters for fit or qa: " << _paramFile << "\n\n"
                << "startHypo: " << _startHypo << "\n\n"
                << "mode: " << _mode << "\n\n"
		<< "number of threads: " << _noOfThreads  << "\n\n"
		<< "ratioMcToData: " << _ratioMcToData  << "\n\n"
		<< "cache amplitudes: " << _cacheAmps  << "\n\n"
		<< "use event weight: " << _useEvtWeight  << "\n\n"
		<< "use phase space hyp: " << _usePhaseSpaceHyp  << "\n\n" 
		<< "pdg table: " << _pdgTableFile << "\n\n" 
            << endl;


      std::vector<std::string>::const_iterator it;

      for (it = _cloneParticle.begin(); it!=_cloneParticle.end(); ++it){
	  std::cout << "clone particles: " << (*it) << "\n";
      }

      for (it=_enabledHyps.begin(); it!=_enabledHyps.end();++it){
          std::cout << "hypothesis\t" << (*it) << "\t enabled\n";
      }
      std::cout << std::endl;
    

      for (it=_mnParFixs.begin(); it!=_mnParFixs.end();++it){
          std::cout << "minuit parameter\t" << (*it) << "\t fixed\n";
      }
      std::cout << std::endl;

      std::cout << "the final state particles are:" << std::endl;
      //      std::vector<std::string>::const_iterator it;
      for (it=_finalStateParticles.begin(); it!=_finalStateParticles.end();++it){
	std::cout << (*it) << "\n";
      }

      std::cout << "\ndecay system:" << std::endl;
      for (it=_decaySystem.begin(); it!=_decaySystem.end();++it){
	std::cout << (*it) << "\n";
      }
      
      std::cout << "\ndecay dynamics:" << std::endl;
      for (it=_dynamics.begin(); it!=_dynamics.end();++it){
	std::cout << (*it) << "\n";
      }
      
      std::cout << "\nreplaced suffix for fit parameter name" << std::endl;
      for (it=_replaceParSuffix.begin(); it!=_replaceParSuffix.end();++it){
	std::cout << (*it) << "\n";
      }

      std::cout << "\nreplaced mass key" << std::endl;
      for (it=_replaceMassKey.begin(); it!=_replaceMassKey.end();++it){
	std::cout << (*it) << "\n";
      }

      std::cout << "\nproduction system:" << std::endl;
      for (it=_productionSystem.begin(); it!=_productionSystem.end();++it){
	std::cout << (*it) << "\n";
      }

      std::cout << "\nproduction formalism:\t" << _productionFormalism << std::endl;


      std::cout << "\nhistograms inv mass for systems" << std::endl;
      for (it=_histMass.begin(); it!=_histMass.end();++it){
	std::cout << (*it) << "\n";
      }
      
      std::cout << "\nhistograms decay angles for systems" << std::endl;
      for (it=_histAngles.begin(); it!=_histAngles.end();++it){
	std::cout << (*it) << "\n";
      }      
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

