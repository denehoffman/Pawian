//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universit??t Bochum 				  //
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

// GlobalEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include <string>
#include <memory>
#include <vector>
#include <map>
#include "AbsChannelEnv.hh"

class ParserBase;
class ParticleTable;
class AbsPawianParameters;

typedef std::vector<std::pair<std::shared_ptr<AbsChannelEnv>, short> > ChannelEnvList;

class GlobalEnv
{

public:
   static GlobalEnv* instance();
   ~GlobalEnv();

  void setup();
   void setup(ParserBase* theParser);
  void replaceParser(ParserBase* theParser);
   void AddEnv(std::shared_ptr<AbsChannelEnv> newEnv, short channelType);
  void setupChannelEnvs();
  void setupChannelEnvs(std::vector<int> channelIDs);
   
  ParticleTable* particleTable() {return _particleTable;}
   ParserBase* parser() {return _theParser;}
   const std::string outputFileNameSuffix() const {return _outputFileNameSuffix;}
   const std::string serializationFileName() const {return _serializationFileName;}
   const std::shared_ptr<AbsChannelEnv> Channel(int id=0) const;
   const std::shared_ptr<AbsChannelEnv> PbarpChannel(int id=0) const;
   const std::shared_ptr<AbsChannelEnv> EpemChannel(int id=0) const;
   const std::shared_ptr<AbsChannelEnv> GammapChannel(int id=0) const;
  const std::shared_ptr<AbsChannelEnv> GGChannel(int id=0) const;
   const std::shared_ptr<AbsChannelEnv> ResChannel(int id=0) const;
  const std::shared_ptr<AbsChannelEnv> PiPiScatteringChannel(int id) const;
   const ChannelEnvList ChannelEnvs() const {return _channelEnvs;}
   const short NoChannels() const {return _channelEnvs.size();}
   std::shared_ptr<AbsPawianParameters> defaultPawianParams();
   std::shared_ptr<AbsPawianParameters> randomPawianParams();
   std::shared_ptr<AbsPawianParameters> startPawianParams() {return _startParams;}
   void setStartPawianParams(std::shared_ptr<AbsPawianParameters> startParams);
   std::vector<std::string> fixedParams();
  std::vector<std::string> paramDependencies();

  std::map<std::string, std::string>& toBeReplacedSuffixMap() {return _toBeReplacedSuffixMap;}
  std::map<std::string, std::string>& alreadyReplacedSuffixMap() {return _alreadyReplacedSuffixMap;}
  std::map<std::string, std::string>& fitParamReplacementMap() {return _fitParamReplacementMap;}
  void fillReplacedSuffixMap(ParserBase* theParser);
  void addIntoToBeReplacedSuffixMap(std::string& toBeReplaced, std::string& replacedBy);
  void addToStringStringMap(const std::string& firstString, const std::string& secondString, std::map<std::string, std::string>& theMap);  
  bool areSuffixMapsIdentical();
  void printFitParameterReplacements();
  std::string topDirPath() const {return _topDirPath;}
  std::string KMatrixPath() const {return _KMatPath;}  
private:
  GlobalEnv();
   static GlobalEnv* _instance;
   bool _alreadySetUp;
   bool _channelEnvsAlredySetup;
   ParticleTable* _particleTable;
   ParserBase* _theParser;
   std::string _outputFileNameSuffix;
   std::string _serializationFileName;
  std::string _topDirPath;
  std::string _KMatPath;

   
  ChannelEnvList _channelEnvs;
  std::shared_ptr<AbsPawianParameters> _startParams;
  std::map<std::string, std::string> _toBeReplacedSuffixMap;
  std::map<std::string, std::string> _alreadyReplacedSuffixMap;
  std::map<std::string, std::string> _fitParamReplacementMap;
};
