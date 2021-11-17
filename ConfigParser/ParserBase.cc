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

// Bertram Kopf (RUB)

#include "ConfigParser/ParserBase.hh"
#include "ErrLogger/ErrLogger.hh"
#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;


ParserBase::ParserBase(int argc,char **argv)
  : _configFile("")
  , _errLogMode(debug)
  , _dataFile("")
  , _mcFile("")
  , _truthFile("")
  , _unitInFile("GEV")
  , _orderInFile("Px Py Pz E")
  , _paramFile("")
  , _serializationFile("")
  , _startHypo("base")
  , _mode("plotmode")
  , _outputFileNameSuffix("")
  , _serverAddress("localhost")
  , _clientNumberWeights("")
  , _verbose(true)
  , _noOfThreads(16)
  , _noOfClients(1)
  , _serverPort(22222)
  , _ratioMcToData(100)
  , _ratioTruthToMc(100)
  , _evoIterations(100)
  , _evoPopulation(20)
  , _evoRatioOfModParams(1.)
  , _cacheAmps(true)
  , _calcContributionError(false)
  , _saveContributionHistos(false)
  , _strErrLogMode("debug")
  , _desc(0)
  , _common(new po::options_description("Common Options"))
  , _config(new po::options_description("Configuration file options"))
  ,_useDataEvtWeight(false)
  ,_useMCEvtWeight(false)
  ,_useTruthEvtWeight(false)
  ,_usePhaseSpaceHyp(false)
  ,_useCohPhaseSpaceHyp(false)
  ,_doScaling(false)
  ,_pdgTableFile("/Particle/pdtNew.table")
  ,_productionDynamics("non")
  ,_productionFormalism("Cano")
  ,_fitqRProduction(false)
  ,_randomSeed(44123)
  ,_genWithModel(true)
  ,_noOfGenEvts(10000)
  ,_noOfDataEvts(1000000)
  ,_noOfBootstrapEvts(10000)
  ,_noOfBootstrapFiles(1)
  ,_boostrapFileName("bootstrapFile")
  ,_tolerance(0.1) 
  ,_noInterScattPoints(0)
  ,_intitial4VecStr("")
  ,_projectile4VecStr("")
  ,_prePathDataFiles("")
  ,_prePathKMatrixFiles("")
  ,_scalingWithChannelID(false)
  ,_fixAllPhases(false)
  ,_nllScalingFactor(1.)
  ,_addChannelScalingId(false)
  ,_channelScalingSuffix("")
  ,_minuitStrategyLevel(1)
  ,_singleChannelId(0)
  ,_sProdExp(0.0)
  ,_useParticleNameForPVecBg(false)
  ,_fixPrimaryChannelScaling(true)
{
  string globalCofigFilePath="/ConfigParser/global.cfg";
  _configFile=getenv("TOP_DIR")+globalCofigFilePath;

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
    ("coupledChannelConfigFile,C",po::value< vector<string> >(&_coupledChannelCfgs),  "Configuration files for coupled channels")
    ("pbarpFiles",po::value< vector<string> >(&_pbarpCfgs),  "Configuration files for pbarp channels")
    ("epemFiles",po::value< vector<string> >(&_epemCfgs),  "Configuration files for epem channels")
    ("resFiles",po::value< vector<string> >(&_resCfgs),  "Configuration files for single resonance decay channels")
    ("ggFiles",po::value< vector<string> >(&_ggCfgs),  "Configuration files for single resonance decay channels")
    ("pipiScatteringFiles",po::value< vector<string> >(&_pipiScatteringCfgs),  "Configuration scattering channels")
    ;
  
  _common->add_options()
    ("errLogMode,e", po::value<string>(&_strErrLogMode)->default_value(_strErrLogMode),"choose mode for Error logger.")
    ("datFile",po::value<string>(&_dataFile), "full path of data file")
    ("mcFile",po::value<string>(&_mcFile), "full path of Monte Carlo file")
    ("truthFile",po::value<string>(&_truthFile), "full path of phasespace generated truth file")
    ("unitInFile",po::value<string>(&_unitInFile),"chosen unit in input files")
    ("orderInFile",po::value<string>(&_orderInFile),"chosen order in input files")
    ("paramFile",po::value<string>(&_paramFile), "file with start parameters for fit or QA (full path)")
    ("serializationFile", po::value<string>(&_serializationFile), "serialized pwa i/o file")
    ("serverAddress", po::value<string>(&_serverAddress), "server address for client mode")
    ("clientNumberWeights", po::value<string>(&_clientNumberWeights), "weights to assign client numbers to channels")
    ("startHypo",po::value<string>(&_startHypo), "choose the hyopthesis to start")
    ("mode",po::value<string>(&_mode), "modes are: pwa, dumpDefaultParams, dumpRandomParams, qaMode, plotmode, spinDensity")
    ("noOfThreads",po::value<int>(&_noOfThreads),  "number of threads for multi threaded mode")
    ("noOfClients",po::value<int>(&_noOfClients),  "number of clients/worker nodes for server mode")
    ("serverPort",po::value<int>(&_serverPort),  "port for client/server mode")
    ("ratioMcToData",po::value<short>(&_ratioMcToData),  "number of MC events defined by ratio #MCs/#Data")
    ("ratioTruthToMc",po::value<short>(&_ratioTruthToMc),  "number of truth events defined by ratio #truths/#Mcs; needed only for qaModeEffCorrection")
    ("evoPopulation",po::value<int>(&_evoPopulation),  "iteration population for evo minimizer")
    ("evoIterations",po::value<int>(&_evoIterations),  "number of iterations for evo minimizer")
    ("evoRatioOfModParams",po::value<double>(&_evoRatioOfModParams),  "chosen (avereaged) ratio of fit parameters to be changed for each population (value between 0. and 1.")
    ("cacheAmps",po::value<bool>(&_cacheAmps),  "cache amplitudes")
    ("contributionError",po::value<bool>(&_calcContributionError),  "calculate the wave contribution error")
    ("saveContributionHistos",po::value<bool>(&_saveContributionHistos),  "creates a histogram root-file for each contribution of Option: calcContribution")
    ("useDataEventWeight",po::value<bool>(&_useDataEvtWeight), "enable/disable input for data event weight")
    ("useMCEventWeight",po::value<bool>(&_useMCEvtWeight), "enable/disable input for Monte Carlo event weight")
    ("useTruthEventWeight",po::value<bool>(&_useTruthEvtWeight), "enable/disable input for truth event weight")
    ("usePhaseSpaceHyp",po::value<bool>(&_usePhaseSpaceHyp), "use hypothesis for phase space")
    ("useCohPhaseSpaceHyp",po::value<bool>(&_useCohPhaseSpaceHyp), "use hypothesis for coherent phase space")
    ("doScaling",po::value<bool>(&_doScaling), "enable/disable prefit with free scaling factor")
    ("name",po::value<string>(&_outputFileNameSuffix), "name that is attached to all otuput file names")
    ("pdgTableFile",po::value<string>(&_pdgTableFile), "path of the pdg-table file relative to the top dir")
    ("randomSeed",po::value<int>(&_randomSeed),  "random seed")
    ("prePathDataFiles",po::value<string>(&_prePathDataFiles), "option to set a pre path to the data files")
    ("prePathKMatrixFiles",po::value<string>(&_prePathKMatrixFiles), "option to set a pre path to the K-matrix files")
    ("scalingWithChannelID",po::value<bool>(&_scalingWithChannelID), "scaling with suffix for channel IDs")
    ("fixAllPhases",po::value<bool>(&_fixAllPhases), "option to fix all phases")
    ("singleChannelId",po::value<unsigned int>(&_singleChannelId), "preferred single channel ID e.g. when running QA")
    ("minuitStrategyLevel",po::value<unsigned int>(&_minuitStrategyLevel),  "set strategy level for minuit fit (1 and 2 are supported)")
    ;
  
  _config->add_options()
    ("verbose",po::value<bool>(&_verbose)->default_value(true), "Determines whether additional information should be emitted")
    ("enableHyp",po::value< vector<string> >(&_enabledHyps), "enable hypotheses")
    ("mnParFix",po::value< vector<string> >(&_mnParFixs),  "minuit parameters can be fixed here")
    ("parameterDependency",po::value< vector<string> >(&_parameterDependencies),  "parameter dependencies")
    ("finalStateParticle",po::value< vector<string> >(&_finalStateParticles),  "name of final state particles")
    ("decay",po::value< vector<string> >(&_decaySystem),  "decay: mother and pair of decay particles")
    ("addDynamics",po::value< vector<string> >(&_dynamics), "add dynamics/line shape for resonances")
    ("replaceParamSuffix",po::value< vector<string> >(&_replaceParSuffix),  "replace suffix for fit parameter name")
    ("replaceMassKey",po::value< vector<string> >(&_replaceMassKey),  "replace Key for the fit parameter of the mass")
    ("replaceProdKey",po::value< vector<string> >(&_replaceProdKey),  "replace key for specific fit parameter of the production")
    ("production",po::value< vector<string> >(&_productionSystem),  "pair of produced particles")
    ("productionDynamics",po::value< string >(&_productionDynamics),  "production dynamics")
    ("productionFormalism",po::value< string >(&_productionFormalism),  "used formalism for the production")
    ("fitqRProduction",po::value<bool>(&_fitqRProduction), "enable/disable fir parameter for individual qR prduction values")
    ("cloneParticle",po::value< vector<string> >(&_cloneParticle),  "particles to be cloned")
    ("preFactor",po::value< vector<string> >(&_preFactor),  "set prefactor for amplitude")
    ("histMass",po::value< vector<string> >(&_histMass),  "histograms inv mass for the selected final state paricles")
    ("histAngles",po::value< vector<string> >(&_histAngles),  "histograms decay angles")
    ("massRangeCuts", po::value< vector<string> > (&_massRangeCuts), "multiple mass range cuts; order: min max particle1 particle2 ...; only events within min and max are accepted")
    ("massRangeAntiCuts", po::value< vector<string> > (&_massRangeAntiCuts), "multiple mass range anti cuts; order: min max particle1 particle2 ...; events within min and max are rejected")
    ("genRange", po::value<string> (&_genRange), "Range of W for generator (to be used if W<=sqrt(s), e.g. in gamma gamma or central production reactions)")
    ("phpGenDynamics", po::value< vector<string> > (&_phpGenDynamics), "dynamics for phase space generated events (only BreitWigner supported so far); order: dynType mass0 width0 particle1 particle2 ...") 
    ("histAngles2D",po::value< vector<string> >(&_histAngles2D),  "2D histogram decay angles")
    ("generateWithModel",po::value<bool>(&_genWithModel),  "generate w/ or w/o model")
    ("noOfGenEvents",po::value<int>(&_noOfGenEvts),  "number of generated events")
    ("noOfDataEvents",po::value<int>(&_noOfDataEvts),  "number of data events for PWA and qa")
    ("noOfBootstrapEvts",po::value<int>(&_noOfBootstrapEvts), "number of required bootstrap events")
    ("noOfBootstrapFiles",po::value<int>(&_noOfBootstrapFiles), "number of bootstrap files")
    ("boostrapFileName", po::value<string>(&_boostrapFileName), "name of the output files from the bootsrap method")
    ("calcContribution",po::value< vector<string> >(&_calcContribution),  "Calculate contribution of partial wave")
    ("minimumTolerance", po::value<double>(&_tolerance), "Minimum tolerance")
    ("noOfInterpolatedScatteringPoints", po::value<int>(&_noInterScattPoints), "number of interpolated scattering points")
    ("kMatrixProdSuffix", po::value< vector<string> >(&_kMatrixProdSuffix),"suffix for fit parameter of the P-vectors in the K-matrix approach")
    ("initial4Vec", po::value<string>(&_intitial4VecStr),"initial 4 Vector (order: E, px, py, pz")
    ("projectile4Vec", po::value<string>(&_projectile4VecStr),"4 Vector of the projectile (order: E, px, py, pz")
    ("nllScalingFactor", po::value<double>(&_nllScalingFactor), "NLL scaling factor for individual channels")
    ("addChannelScalingId",po::value<bool>(&_addChannelScalingId), "add Id to parameter name for channel scaling")
    ("channelScalingSuffix",po::value<string>(& _channelScalingSuffix), "suffix for name of channel scling parameter")
    ("sProdExponent", po::value<double>(&_sProdExp)->default_value(_sProdExp),"exponent for the s dependency in the production (supported for FVectorIntensity and Gamgam fits")
    ("useParticleNameForPVecBg",po::value<bool>(&_useParticleNameForPVecBg), "yes: (cloned) particle name; no: K-matrix name")
    ("fixPrimaryChannelScaling",po::value<bool>(&_fixPrimaryChannelScaling), "fix/release scaling of first channel for coupled channel mode")
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
	strError << "Error accessing configuration file " << _configFile;
	std::cout << cmdline_options << endl;

	throw runtime_error(strError.str());
      }


    store(po::parse_config_file(ifs, config_file_options), vm);
    po::notify(vm);

    // Emit a help message, if necessary
    if (vm.count("help"))
    {
      std::cout << config_file_options << endl;
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
      WarningMsg << "ErrorLogger not (properly) set -> Use mode 'DEBUG' " ;  // << endmsg;
    }

    if (_evoRatioOfModParams<=0. || _evoRatioOfModParams>1.){
      Alert << "_evoRatioOfModParams = " << _evoRatioOfModParams << " not possible\n"
    	    << "value must be between 0. and 1. !!!!" << endmsg;
      exit(1);
    }

    if(_verbose){
      std::cout << "\nRunning with the following options using " << _configFile << ":\n\n"
                << "Error log mode: " << _errLogMode <<"\n\n"
                << "data file: " << _dataFile <<"\n\n"
                << "mc file: " << _mcFile <<"\n\n"
		<< "truth file: " << _truthFile <<"\n\n"
		<< "chosen unit in input file: " << _unitInFile << "\n\n"
		<< "chosen order in file: " << _orderInFile << "\n\n"
                << "file with start parameters for fit or qa: " << _paramFile << "\n\n"
                << "startHypo: " << _startHypo << "\n\n"
                << "mode: " << _mode << "\n\n"
		<< "number of threads: " << _noOfThreads  << "\n\n"
		<< "ratioMcToData: " << _ratioMcToData  << "\n\n"
		<< "ratioTruthToMc: " << _ratioTruthToMc  << "\n\n"
		<< "cache amplitudes: " << _cacheAmps  << "\n\n"
		<< "use data event weight: " << _useDataEvtWeight  << "\n\n"
		<< "use Monte Carlo event weight: " << _useMCEvtWeight  << "\n\n"
		<< "use truth event weight: " << _useTruthEvtWeight  << "\n\n"
		<< "use phase space hyp: " << _usePhaseSpaceHyp  << "\n\n"
                << "use coherent phase space hyp: " << _useCohPhaseSpaceHyp  << "\n\n"
		<< "prefit with free scaling factor: " << _doScaling << "\n\n"
		<< "pdg table: " << _pdgTableFile << "\n\n"
	        << "minimumTolerance: " << _tolerance << "\n\n"
                << "no of interpolated scattering points: " << _noInterScattPoints << "\n\n"
		<< "initial 4vector: " << _intitial4VecStr << "\n\n"
		<< "projectile 4vector: " << _projectile4VecStr << "\n\n"
		<< "prePathDataFiles: " << _prePathDataFiles << "\n\n"
		<< "prePathKMatrixFiles: " << _prePathKMatrixFiles << "\n\n"
		<< "scalingWithChannelID: " <<  _scalingWithChannelID << "\n\n"
		<< "fixAllPhases: " << _fixAllPhases << "\n\n"
    << "singleChannelId: " << _singleChannelId << "\n\n"
		<< "nllScalingFactor: " << _nllScalingFactor << "\n\n"
		<< "addChannelScalingId: " <<_addChannelScalingId << "\n\n"
		<< "channelScalingSuffix: " <<_channelScalingSuffix << "\n\n"
		<< "minuitStrategyLevel: " << _minuitStrategyLevel << "\n\n"
	        << "useParticleNameForPVecBg: " << _useParticleNameForPVecBg << "\n\n"
                << endl;


      std::vector<std::string>::const_iterator it;

      for (it = _coupledChannelCfgs.begin(); it!=_coupledChannelCfgs.end(); ++it){
	  std::cout << "Coupled channel configuration file: " << (*it) << "\n";
      }
      for (it = _pbarpCfgs.begin(); it!=_pbarpCfgs.end(); ++it){
         std::cout << "pbarp channel configuration file: " << (*it) << "\n";
      }
      for (it = _epemCfgs.begin(); it!=_epemCfgs.end(); ++it){
         std::cout << "epem channel configuration file: " << (*it) << "\n";
      }
      for (it = _resCfgs.begin(); it!=_resCfgs.end(); ++it){
         std::cout << "res channel configuration file: " << (*it) << "\n";
      }
      for (it = _ggCfgs.begin(); it!=_ggCfgs.end(); ++it){
         std::cout << "gg channel configuration file: " << (*it) << "\n";
      }
      for (it = _pipiScatteringCfgs.begin(); it!=_pipiScatteringCfgs.end(); ++it){
         std::cout << "pipi scattering channel configuration file: " << (*it) << "\n";
      }
      for (it = _cloneParticle.begin(); it!=_cloneParticle.end(); ++it){
	  std::cout << "clone particles: " << (*it) << "\n";
      }
      for (it = _preFactor.begin(); it!=_preFactor.end(); ++it){
	  std::cout << "preFactors: " << (*it) << "\n";
      }
      for (it=_enabledHyps.begin(); it!=_enabledHyps.end();++it){
          std::cout << "hypothesis\t" << (*it) << "\t enabled\n";
      }
      std::cout << std::endl;


      for (it=_mnParFixs.begin(); it!=_mnParFixs.end();++it){
          std::cout << "minuit parameter\t" << (*it) << "\t fixed\n";
      }
      std::cout << std::endl;
      for (it=_parameterDependencies.begin(); it!=_parameterDependencies.end();++it){
	std::cout << "parameter dependencies:\t" << (*it) << "\n";
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

      std::cout << "\nreplaced prod key" << std::endl;
      for (it=_replaceProdKey.begin(); it!=_replaceProdKey.end();++it){
	std::cout << (*it) << "\n";
      }

      std::cout << "\nproduction system:" << std::endl;
      for (it=_productionSystem.begin(); it!=_productionSystem.end();++it){
	std::cout << (*it) << "\n";
      }

      std::cout << "\nK-matrix prod suffix:" << std::endl;
      for (it=_kMatrixProdSuffix.begin(); it!=_kMatrixProdSuffix.end();++it){
	std::cout << (*it) << "\n";
      }

      std::cout << "\nproduction formalism:\t" << _productionFormalism << std::endl;
      std::cout << "fitqRProduction:\t" << _fitqRProduction << std::endl;


      std::cout << "\nhistograms inv mass for systems" << std::endl;
      for (it=_histMass.begin(); it!=_histMass.end();++it){
	std::cout << (*it) << "\n";
      }

      std::cout << "\nhistograms decay angles for systems" << std::endl;
      for (it=_histAngles.begin(); it!=_histAngles.end();++it){
	std::cout << (*it) << "\n";
      }

      // std::cout << "\nmass range" << std::endl;
      // std::cout << _massRange << "\n";

      std::cout << "\nmass range cuts:" << std::endl;
      for (it = _massRangeCuts.begin(); it!=_massRangeCuts.end(); ++it){
	  std::cout << (*it) << "\n";
      }

      std::cout << "\nmass range anti-cuts:" << std::endl;
      for (it = _massRangeAntiCuts.begin(); it!=_massRangeAntiCuts.end(); ++it){
	  std::cout << (*it) << "\n";
      }

      std::cout << "\ngenerator range:\t" << _genRange << std::endl;

      std::cout << "\nphp dynamics:" << std::endl;
      for (it = _phpGenDynamics.begin(); it!=_phpGenDynamics.end(); ++it){
	  std::cout << (*it) << "\n";
      }

      std::cout << "\n2Dhistogram decay angles for systems" << std::endl;
      for (it=_histAngles2D.begin(); it!=_histAngles2D.end();++it){
	std::cout << (*it) << "\n";
      }
      std::cout << "\nrandom seed:\t" << _randomSeed << std::endl;
      std::cout << "\ngenerate with model:\t" << _genWithModel << std::endl;
      std::cout << "\nnumber of generated events:\t" << _noOfGenEvts << std::endl;
      std::cout << "\nnumber of data events:\t" << _noOfDataEvts << std::endl;
      std::cout << "\nnumber of bootstrap events:\t" << _noOfBootstrapEvts << std::endl;
      std::cout << "\nnumber of bootstrap files:\t" << _noOfBootstrapFiles << std::endl;
      std::cout << "\nboostrap file name:\t" << _boostrapFileName << std::endl;
      std::cout << "sProdExponent =" << _sProdExp << std::endl;
    }

  }


  catch( std::exception & e )
  {
    cerr << "Error parsing the command line:" << endl;
    cerr << e.what() << std::endl;
    cerr << "You can use -h or --help to obtain the description of the program parameters." << endl;
    //    cerr << "This is the command line options\n" << endl;
    exit(1);
    //    return false;
  }
  catch(...){
      std::cerr << "Error parsing the command line. Use -h or --help to see the description of the program paramters." << endl;
      return false;
  }

  return true;
}

