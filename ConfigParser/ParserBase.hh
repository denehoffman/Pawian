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
  const std::string truthFile() const {return _truthFile;}
  const std::string unitInFile() const {return _unitInFile;}
  const std::string orderInFile() const {return _orderInFile;}
  const std::string fitParamFile() const {return _paramFile;}
  const std::string serializationFile() const{ return _serializationFile; }
  const std::string outputFileNameSuffix() const {return _outputFileNameSuffix;}
  const std::string serverAddress() const {return _serverAddress;}
  const std::string clientNumberWeights() const {return _clientNumberWeights;}
  const std::vector<std::string>& enabledHyps() const { return _enabledHyps; }
  const std::string startHypo() const {return _startHypo;}
  const std::string mode() const {return _mode;}
  const std::vector<std::string>& fixedParams() const { return _mnParFixs; }
  const std::vector<std::string>& parameterDependencies() const { return _parameterDependencies; }
  const int noOfThreads() const {return _noOfThreads;}
  const int noOfClients() const {return _noOfClients;}
  const int serverPort() const {return _serverPort;}
  const short ratioMcToData() const {return _ratioMcToData;}
  const short ratioTruthToMc() const {return _ratioTruthToMc;}
  const int evoPopulation() const {return _evoPopulation;}
  const int evoIterations() const {return _evoIterations;}
  const double evoRatioOfModParams() const {return _evoRatioOfModParams;}
  const bool cacheAmps() const {return _cacheAmps;}
  const bool calcContributionError() const {return _calcContributionError;}
  const bool saveContributionHistos() const {return _saveContributionHistos;}
  const bool useDataEvtWeight() const {return _useDataEvtWeight; }
  const bool useMCEvtWeight() const {return _useMCEvtWeight; }
  const bool useTruthEvtWeight() const {return _useTruthEvtWeight; }
  const bool usePhaseSpaceHyp() const {return _usePhaseSpaceHyp; }
  const bool useCohPhaseSpaceHyp() const {return _useCohPhaseSpaceHyp; }
  const bool doScaling() const {return _doScaling;}
  const std::string pdgTableFile() const {return _pdgTableFile;}
  const std::vector<std::string>& coupledChannelCfgs() const { return _coupledChannelCfgs; }
  const std::vector<std::string>& pbarpCfgs() const {return _pbarpCfgs; }
  const std::vector<std::string>& epemCfgs() const {return _epemCfgs; }
  const std::vector<std::string>& resCfgs() const {return _resCfgs; }
  const std::vector<std::string>& ggCfgs() const {return _ggCfgs; }
  const std::vector<std::string>& pipiScatteringCfgs() const {return _pipiScatteringCfgs; }
  const std::vector<std::string>& finalStateParticles() const { return _finalStateParticles; }
  const std::vector<std::string>& decaySystem() const { return _decaySystem; }
  const std::vector<std::string>& replaceSuffixNames() const { return _replaceParSuffix; }
  const std::vector<std::string>& replaceMassKey() const { return _replaceMassKey; }
  const std::vector<std::string>& replaceProdKey() const { return _replaceProdKey; }
  const std::vector<std::string>& decayDynamics() const { return _dynamics;}
  const std::vector<std::string>& productionSystem() const { return _productionSystem; }
  const std::string productionDynamics() const {return _productionDynamics;}
  const std::string productionFormalism() const { return _productionFormalism; }
   const bool fitqRProduction() const {return _fitqRProduction;}
  const std::vector<std::string>& cloneParticle() const { return _cloneParticle; }
  const std::vector<std::string>& preFactor() const { return _preFactor;}
  const std::vector<std::string>& histMassNames() const { return _histMass; }
  const std::vector<std::string>& histAngleNames() const { return _histAngles;}
  const std::vector<std::string>& histAngleNames2D() const { return _histAngles2D;}
  const std::vector<std::string>& calcContribution() const { return _calcContribution;}
  const std::vector<std::string>& massRangeCuts() const {return _massRangeCuts;}
  const std::vector<std::string>& massRangeAntiCuts() const {return _massRangeAntiCuts;}
  const std::string& genRange() const {return _genRange;}
  const std::vector<std::string>& phpGenDynamics() const {return _phpGenDynamics;}
  const int randomSeed() const {return _randomSeed;}
  const bool generateWithModel() const {return _genWithModel;}
  const int noOfGenEvts() const {return _noOfGenEvts;}
  const int noOfDataEvts() const {return _noOfDataEvts;}
  const int noOfBootstrapEvts() const {return _noOfBootstrapEvts;}
  const int noOfBootstrapFiles() const {return _noOfBootstrapFiles;}
  const std::string bootstrapFileName() const {return _boostrapFileName;}
  const double tolerance() const {return _tolerance;}
  const int noOfInterpolatedScattPoints() const {return _noInterScattPoints;}
  const std::vector<std::string>& addKmatrixProdSuffix() const {return _kMatrixProdSuffix;}
  const std::string intitial4Vec() const {return _intitial4VecStr;}
  const std::string projectile4Vec() const {return _projectile4VecStr;}
  const std::string prePathDataFiles() const {return _prePathDataFiles;}
  const std::string prePathKMatrixFiles() const {return _prePathKMatrixFiles;}
  const bool scalingWChannelID() const {return _scalingWithChannelID;}
  const bool fixAllPhases() const {return _fixAllPhases;}
  const double nllScalingFactor() const {return _nllScalingFactor;}
  const bool addChannelScalingId() const{return _addChannelScalingId;}
  const std::string channelScalingSuffix() const {return _channelScalingSuffix;}
  const unsigned int minuitStrategyLevel() const{return _minuitStrategyLevel;}
  const unsigned int singleChannelId() const{ return _singleChannelId;}
  const double sProdExponent() const {return _sProdExp;}
  const bool useParticleNameForPVecBg() const {return _useParticleNameForPVecBg;}
  const bool fixPrimaryChannelScaling() const {return _fixPrimaryChannelScaling;}
  const unsigned int prodOrbMom1DFit() const {return _prodOrbMom1DFit;}
protected:
  virtual bool parseCommandLine(int argc,char **argv); 
  
  std::string _configFile;
  enErrLogMode _errLogMode;
  std::string _dataFile;
  std::string _mcFile;
  std::string _truthFile;
  std::string _unitInFile;
  std::string _orderInFile;
  std::string _paramFile;
  std::string _serializationFile;
  std::string _startHypo;
  std::string _mode;		  
  std::string _outputFileNameSuffix;
  std::string _serverAddress;
  std::string _clientNumberWeights;
  std::vector<std::string> _enabledHyps;
  std::vector<std::string> _mnParFixs;
  std::vector<std::string> _parameterDependencies;
  bool _verbose;
  int _noOfThreads;
  int _noOfClients;
  int _serverPort;
  short _ratioMcToData;
  short _ratioTruthToMc;
  int _evoIterations;
  int _evoPopulation;
  double _evoRatioOfModParams;
  bool _cacheAmps;
  bool _calcContributionError;
  bool _saveContributionHistos;
  std::string _strErrLogMode;
  po::options_description* _desc;
  po::options_description* _common;
  po::options_description* _config;
  bool _useDataEvtWeight;
  bool _useMCEvtWeight;
  bool _useTruthEvtWeight;
  bool _usePhaseSpaceHyp;
  bool _useCohPhaseSpaceHyp;
  bool _doScaling;
  std::string _pdgTableFile;
  std::vector<std::string> _coupledChannelCfgs;
  std::vector<std::string> _pbarpCfgs;
  std::vector<std::string> _epemCfgs;
  std::vector<std::string> _resCfgs;
  std::vector<std::string> _ggCfgs;
  std::vector<std::string> _pipiScatteringCfgs;
  std::vector<std::string> _finalStateParticles;
  std::vector<std::string> _decaySystem;
  std::vector<std::string> _dynamics;
  std::vector<std::string> _replaceParSuffix;
  std::vector<std::string>  _replaceMassKey;
  std::vector<std::string>  _replaceProdKey;
  std::vector<std::string> _productionSystem;
  std::string _productionDynamics;
  std::string _productionFormalism;
  bool _fitqRProduction;
  std::vector<std::string> _cloneParticle;
  std::vector<std::string> _preFactor;
  std::vector<std::string> _histMass;
  std::vector<std::string> _histAngles;
  std::vector<std::string>  _massRangeCuts;
  std::vector<std::string>  _massRangeAntiCuts;
  std::string  _genRange;
  std::vector<std::string>  _phpGenDynamics;
  std::vector<std::string> _histAngles2D;
  std::vector<std::string> _calcContribution;
  int _randomSeed;
  bool _genWithModel;
  int _noOfGenEvts;
  int _noOfDataEvts;
  int _noOfBootstrapEvts;
  int _noOfBootstrapFiles;
  std::string _boostrapFileName;
  double _tolerance;
  int _noInterScattPoints;
  std::vector<std::string> _kMatrixProdSuffix;
  std::string _intitial4VecStr;
  std::string _projectile4VecStr;
  std::string _prePathDataFiles;
  std::string _prePathKMatrixFiles;
  bool _scalingWithChannelID;
  bool _fixAllPhases;
  double _nllScalingFactor;
  bool _addChannelScalingId;
  std::string _channelScalingSuffix;
  unsigned int _minuitStrategyLevel;
  unsigned int _singleChannelId;
  double _sProdExp;
  bool _useParticleNameForPVecBg;
  bool _fixPrimaryChannelScaling;
  unsigned int _prodOrbMom1DFit;
};
