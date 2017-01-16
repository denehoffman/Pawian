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

#include "PwaUtils/GlobalEnv.hh"
#include "ConfigParser/ParserBase.hh"
#include "PwaUtils/AbsLh.hh"
#include "Particle/PdtParser.hh"
#include "ErrLogger/ErrLogger.hh"

#include "FitParams/AbsPawianParameters.hh"
#include "FitParams/ParamFactory.hh"


GlobalEnv* GlobalEnv::_instance=0;

GlobalEnv* GlobalEnv::instance(){
   if(0==_instance) _instance = new GlobalEnv();
   return _instance;
}


GlobalEnv::GlobalEnv() :
   _alreadySetUp(false) ,
   _channelEnvsAlredySetup(false),
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
}

void GlobalEnv::setupChannelEnvs(){
  if(_channelEnvsAlredySetup){
    Alert << "channel environments already setup!!!" << endmsg;
    exit(1);
  }

  int id=0;
  for(auto it = _channelEnvs.begin(); it!=_channelEnvs.end();++it){
      (*it).first->setup(id);
      if (!(*it).first->checkReactionChain()){
	Alert << "Something wrong with the reaction chain for channelTypeName: " << (*it).first->channelTypeName() << endmsg;
	exit(1);
      }
      (*it).first->setWignerDRefs();
      ++id;
   }
  _channelEnvsAlredySetup=true;
}

std::vector<std::string> GlobalEnv::fixedParams(){
  std::vector<std::string> result;
  std::vector<std::string>::iterator itStr;
  for(auto it = _channelEnvs.begin(); it!=_channelEnvs.end();++it){
   std::vector<std::string> currentFixedParmaList=(*it).first->parser()->fixedParams();
   for(itStr=currentFixedParmaList.begin(); itStr!=currentFixedParmaList.end(); ++itStr){
 //   std::vector<std::string>::iterator itFindStr = result.find(*itStr);
 //    if (itFindStr == result.end()) result.push_back(*itStr);
   result.push_back(*itStr);
  }
 }
 return result;
}

std::shared_ptr<AbsPawianParameters> GlobalEnv::defaultPawianParams(){
  std::shared_ptr<AbsPawianParameters> result=ParamFactory::instance()->getParametersPointer("Pawian");
  for(auto it = _channelEnvs.begin(); it!=_channelEnvs.end();++it){
    if(0==(*it).first->Lh()){
      Alert << "for channel id: " << (*it).first->channelID() << " typeName: " << (*it).first->channelTypeName() << endmsg; 
      Alert << "likelihood object is a NULL pou=inter!!!!" << endmsg;
      exit(1);
    }
    (*it).first->Lh()->fillDefaultParams(result);
   }
  return result;
}

void GlobalEnv::setStartPawianParams(std::shared_ptr<AbsPawianParameters> startParams){
  std::shared_ptr<AbsPawianParameters> defaultParams=GlobalEnv::instance()->defaultPawianParams();
  _startParams=defaultParams->paramsWithSameOrder(startParams);
}

void GlobalEnv::setup(ParserBase* theParser){
   if(_alreadySetUp){
      Alert << "GlobalEnv already set up!" << endmsg;
      exit(1);
   }

   _alreadySetUp = true;

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

void GlobalEnv::replaceParser(ParserBase* theParser){
  _alreadySetUp = false;
  InfoMsg << "Now replace the parser!!!" << endmsg;
  setup(theParser);  
}

void GlobalEnv::addIntoToBeReplacedSuffixMap(std::string& toBeReplaced, std::string& replacedBy){

  //check if replacement still does not exist
  std::map<std::string, std::string>::iterator it;
  for(it=_toBeReplacedSuffixMap.begin(); it!=_toBeReplacedSuffixMap.end(); ++it){
    if(it->first == toBeReplaced){
      if(it->second != replacedBy){
	Alert <<"ambiguous replacement of replaceParamSuffixes!!!" 
	      << it->first << " cannot replaced by " << replacedBy 
	      <<"/nsince replacement of " <<  it->first << " by " << it->second << " already exists!!!" << endmsg;
	exit(1); 
      }
      else return;
    }
  }

  _toBeReplacedSuffixMap[toBeReplaced]=replacedBy; 
}

bool GlobalEnv::areSuffixMapsIdentical(){
  //check if _toBeReplacedSuffixMap ==  _alreadyReplacedSuffixMap
  std::map<std::string, std::string>::iterator itToFind;
  std::map<std::string, std::string>::iterator it;
  for(it=_toBeReplacedSuffixMap.begin(); it!=_toBeReplacedSuffixMap.end(); ++it){
    itToFind = _alreadyReplacedSuffixMap.find(it->first);
    if(itToFind != _alreadyReplacedSuffixMap.end()){
      DebugMsg << "found replacement of parameter suffix from " << it->first << " to " << it->second << endmsg; 
    }
    else{
      Alert << "It is not possible to replace of parameter suffix " << it->first << " by " << it->second << endmsg;
      Alert << "the already replaced map contains " << endmsg;
      for(itToFind=_alreadyReplacedSuffixMap.begin(); itToFind!=_alreadyReplacedSuffixMap.end(); ++itToFind){
	Alert << itToFind->first << "\t" << itToFind->second << endmsg;
      }

      std::map<std::string, std::string>::iterator it2;
     Alert << "\nthe to be replaced map contains " << endmsg;
     for(it2=_toBeReplacedSuffixMap.begin(); it2!=_toBeReplacedSuffixMap.end(); ++it2){
       Alert << it2->first << "\t" << it2->second << endmsg;
     }
      exit(1);
    }
  }
  return true;
}

void GlobalEnv::addToStringStringMap(const std::string& firstString, const std::string& secondString, std::map<std::string, std::string>& theMap){
  //check if replacement still does not exist
  std::map<std::string, std::string>::iterator it;
  for(it=theMap.begin(); it!=theMap.end(); ++it){
    if(it->first == firstString) return;
  }
  theMap[firstString]=secondString;
}

void GlobalEnv::printFitParameterReplacements(){
  InfoMsg << "\n*** The fit parameter replacements are ***" << endmsg;
  std::map<std::string, std::string>::iterator it;
  for(it=_fitParamReplacementMap.begin(); it!=_fitParamReplacementMap.end(); ++it){
    InfoMsg << it->first << "\treplaced by\t" << it->second << endmsg;
  }
}  
