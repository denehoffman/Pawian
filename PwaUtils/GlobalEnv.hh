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

// GlobalEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include <string>
#include <memory>
#include "AbsChannelEnv.hh"

class ParserBase;
class ParticleTable;
class FitParColBase;

typedef std::vector<std::pair<std::shared_ptr<AbsChannelEnv>, short> > ChannelEnvList;

class GlobalEnv
{

public:
   static GlobalEnv* instance();
   GlobalEnv();
   ~GlobalEnv();

   void setup(ParserBase* theParser);
  void replaceParser(ParserBase* theParser);
   void AddEnv(std::shared_ptr<AbsChannelEnv> newEnv, short channelType);

   ParticleTable* particleTable() {return _particleTable;}
   ParserBase* parser() {return _theParser;}
   std::shared_ptr<FitParColBase> fitParColBase(){ return _theFitParColBase;}
   const std::string outputFileNameSuffix() const {return _outputFileNameSuffix;}
   const std::string serializationFileName() const {return _serializationFileName;}
   const std::shared_ptr<AbsChannelEnv> Channel(int id=0) const;
   const std::shared_ptr<AbsChannelEnv> PbarpChannel(int id=0) const;
   const std::shared_ptr<AbsChannelEnv> EpemChannel(int id=0) const;
   const std::shared_ptr<AbsChannelEnv> GammapChannel(int id=0) const;
   const std::shared_ptr<AbsChannelEnv> ResChannel(int id=0) const;
   const ChannelEnvList ChannelEnvs() const {return _channelEnvs;}
   void CreateDefaultParameterSet();
   const short NoChannels() const {return _channelEnvs.size();}
   fitParCol DefaultParamVal(){ return _defaultVal;}
   fitParCol DefaultParamErr(){ return _defaultErr;}

private:
   static GlobalEnv* _instance;
   bool _alreadySetUp;
   ParticleTable* _particleTable;
   ParserBase* _theParser;
   std::shared_ptr<FitParColBase> _theFitParColBase;
   std::string _outputFileNameSuffix;
   std::string _serializationFileName;

   ChannelEnvList _channelEnvs;

   fitParCol _defaultVal;
   fitParCol _defaultErr;
};
