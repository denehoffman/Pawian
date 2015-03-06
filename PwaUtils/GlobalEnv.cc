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

#include "PwaUtils/GlobalEnv.hh"
#include "ConfigParser/ParserBase.hh"
#include "PwaUtils/AbsLh.hh"
#include "Particle/PdtParser.hh"
#include "ErrLogger/ErrLogger.hh"


GlobalEnv* GlobalEnv::_instance=0;

GlobalEnv* GlobalEnv::instance(){
   if(0==_instance) _instance = new GlobalEnv();
   return _instance;
}


GlobalEnv::GlobalEnv() :
   _alreadySetUp(false) ,
   _theParser(0)
{
}



const std::shared_ptr<AbsChannelEnv> GlobalEnv::Channel(int id) const
{
   return _channelEnvs.at(id).first;
}

const std::shared_ptr<AbsChannelEnv> GlobalEnv::PbarpChannel(int id) const
{

   if(_channelEnvs.at(id).second != AbsChannelEnv::CHANNEL_PBARP){
      Alert << "Faultily accessing non-pbarp channel environment." << endmsg;
   }

   return _channelEnvs.at(id).first;
}



const std::shared_ptr<AbsChannelEnv> GlobalEnv::EpemChannel(int id) const
{

   if(_channelEnvs.at(id).second != AbsChannelEnv::CHANNEL_EPEM){
      Alert << "Faultily accessing non-epem channel environment." << endmsg;
   }

   return _channelEnvs.at(id).first;
}

const std::shared_ptr<AbsChannelEnv> GlobalEnv::GammapChannel(int id) const
{

   if(_channelEnvs.at(id).second != AbsChannelEnv::CHANNEL_GAMMAP){
      Alert << "Faultily accessing non-epem channel environment." << endmsg;
   }

   return _channelEnvs.at(id).first;
}


const std::shared_ptr<AbsChannelEnv> GlobalEnv::ResChannel(int id) const
{

   if(_channelEnvs.at(id).second != AbsChannelEnv::CHANNEL_RES){
      Alert << "Faultily accessing non-res channel environment." << endmsg;
   }

   return _channelEnvs.at(id).first;
}



void GlobalEnv::AddEnv(std::shared_ptr<AbsChannelEnv> newEnv, short envType){
   _channelEnvs.push_back(std::pair<std::shared_ptr<AbsChannelEnv>, short>(newEnv, envType));
   newEnv->setup(_channelEnvs.size() - 1);
}



void GlobalEnv::CreateDefaultParameterSet(){
   fitParCol defaultVal;
   fitParCol defaultErr;
   for(auto it = _channelEnvs.begin(); it!=_channelEnvs.end();++it){
      (*it).first->Lh()->getDefaultParams(defaultVal, defaultErr);
   }
   _defaultVal = defaultVal;
   _defaultErr = defaultErr;
}



void GlobalEnv::setup(ParserBase* theParser){
   if(_alreadySetUp){
      Alert << "GlobalEnv already set up!" << endmsg;
      exit(1);
   }

   _alreadySetUp = true;

   _theFitParColBase = std::shared_ptr<FitParColBase>(new FitParColBase());
   _theParser = theParser;
   _outputFileNameSuffix = theParser->outputFileNameSuffix();
   _serializationFileName = theParser->serializationFile();

   // pdtTable
   PdtParser pdtParser;
   std::string theSourcePath=getenv("TOP_DIR");
   std::string pdtFileRelPath=theParser->pdgTableFile();
   std::string pdtFile(theSourcePath+pdtFileRelPath);
   _particleTable = new ParticleTable;

   if (!pdtParser.parse(pdtFile, *_particleTable)) {
      Alert << "can not parse particle table " << pdtFile << endmsg;
      exit(1);
   }

}
