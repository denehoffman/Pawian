
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

#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

#include <boost/version.hpp>

#if BOOST_VERSION < 103600
#error "Error: Boost should at least have version 1.36 !"
#endif // BOOST_VERSION

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#pragma once


namespace po = boost::program_options;

class ParserBase {
public:
  typedef enum tagerrLogMode { debug, trace, routine, warning, error, alert } enErrLogMode;

  ParserBase(int argc,char **argv);
  virtual ~ParserBase(){;}
    
  const std::string& getConfigFile() const { return _configFile;}
  const enErrLogMode& getErrLogMode() const { return _errLogMode; }
  const std::string dataFile() const {return _dataFile;}
  const std::string mcFile() const {return _mcFile;}
  const std::string unitInFile() const {return _unitInFile;}
  const std::string orderInFile() const {return _orderInFile;}
  const std::string fitParamFile() const {return _paramFile;}
  const std::string serializationFile() const{ return _serializationFile; }
  const std::string outputFileNameSuffix() const {return _outputFileNameSuffix;}
  const std::string serverAddress() const {return _serverAddress;}
  const std::vector<std::string>& enabledHyps() const { return _enabledHyps; }
  const std::string startHypo() const {return _startHypo;}
  const std::string mode() const {return _mode;}
  const std::vector<std::string>& fixedParams() const { return _mnParFixs; } 
  const int noOfThreads() const {return _noOfThreads;}
  const int noOfClients() const {return _noOfClients;}
  const int serverPort() const {return _serverPort;}
  const int ratioMcToData() const {return _ratioMcToData;}
  const int evoPopulation() const {return _evoPopulation;}
  const int evoIterations() const {return _evoIterations;}
  const bool cacheAmps() const {return _cacheAmps;}
  const bool calcContributionError() const {return _calcContributionError;}
  const bool saveContributionHistos() const {return _saveContributionHistos;}
  const bool useEvtWeight() const {return _useEvtWeight; }
  const bool usePhaseSpaceHyp() const {return _usePhaseSpaceHyp; }
  const std::string pdgTableFile() const {return _pdgTableFile;}
  const std::vector<std::string>& coupledChannelCfgs() const { return _coupledChannelCfgs; }
  const std::vector<std::string>& finalStateParticles() const { return _finalStateParticles; }
  const std::vector<std::string>& decaySystem() const { return _decaySystem; }
  const std::vector<std::string>& replaceSuffixNames() const { return _replaceParSuffix; }
  const std::vector<std::string>& replaceMassKey() const { return _replaceMassKey; }
  const std::vector<std::string>& decayDynamics() const { return _dynamics;}
  const std::vector<std::string>& productionSystem() const { return _productionSystem; }
  const std::string productionFormalism() const { return _productionFormalism; }
  const bool useProductionBarrier() const { return _useProductionBarrier; }
  const std::vector<std::string>& cloneParticle() const { return _cloneParticle; }
  const std::vector<std::string>& preFactor() const { return _preFactor;}
  const std::vector<std::string>& histMassNames() const { return _histMass; }
  const std::vector<std::string>& histAngleNames() const { return _histAngles;}
  const std::string& massRange() const { return _massRange;} 
  const std::vector<std::string>& histAngleNames2D() const { return _histAngles2D;}
  const std::vector<std::string>& calcContribution() const { return _calcContribution;}

  const int randomSeed() const {return _randomSeed;}
  const bool generateWithModel() const {return _genWithModel;}
  const int noOfGenEvts() const {return _noOfGenEvts;}
  const int noOfDataEvts() const {return _noOfDataEvts;}
  const double tolerance() const {return _tolerance;}
 
protected:
  virtual bool parseCommandLine(int argc,char **argv); 

  std::string _configFile;
  enErrLogMode _errLogMode;
  std::string _dataFile;
  std::string _mcFile;
  std::string _unitInFile;
  std::string _orderInFile;
  std::string _paramFile;
  std::string _serializationFile;
  std::string _startHypo;
  std::string _mode;		  
  std::string _outputFileNameSuffix;
  std::string _serverAddress;
  std::vector<std::string> _enabledHyps;
  std::vector<std::string> _mnParFixs;
  bool _verbose;
  int _noOfThreads;
  int _noOfClients;
  int _serverPort;
  int _ratioMcToData;
  int _evoIterations;
  int _evoPopulation;
  bool _cacheAmps;
  bool _calcContributionError;
  bool _saveContributionHistos;
  std::string _strErrLogMode;
  po::options_description* _desc;
  po::options_description* _common;
  po::options_description* _config;
  bool _useEvtWeight;
  bool _usePhaseSpaceHyp;
  std::string _pdgTableFile;
  std::vector<std::string> _coupledChannelCfgs;
  std::vector<std::string> _finalStateParticles;
  std::vector<std::string> _decaySystem;
  std::vector<std::string> _dynamics;
  std::vector<std::string> _replaceParSuffix;
  std::vector<std::string>  _replaceMassKey;
  std::vector<std::string> _productionSystem;
  std::string _productionFormalism;
  bool _useProductionBarrier;
  std::vector<std::string> _cloneParticle;
  std::vector<std::string> _preFactor;
  std::vector<std::string> _histMass;
  std::vector<std::string> _histAngles;
  std::string  _massRange;
  std::vector<std::string> _histAngles2D;
  std::vector<std::string> _calcContribution;
  int _randomSeed;
  bool _genWithModel;
  int _noOfGenEvts;
  int _noOfDataEvts;
  double _tolerance;
};
