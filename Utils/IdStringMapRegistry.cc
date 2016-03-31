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

// IdStringMapRegistry class definition file. -*- C++ -*-
// Copyright 2016 Bertram Kopf

#include "Utils/IdStringMapRegistry.hh"
#include "ErrLogger/ErrLogger.hh"


IdStringMapRegistry* IdStringMapRegistry::_instance=0;

IdStringMapRegistry* IdStringMapRegistry::instance()
{
  if (0==_instance) _instance = new IdStringMapRegistry();
  return _instance;
}

IdStringMapRegistry::IdStringMapRegistry() :
  _stringPairMapCounter(0)
  ,_stringMapCounter(0)
{
}

IdStringMapRegistry::~IdStringMapRegistry()
{
}

unsigned short IdStringMapRegistry::stringStringId(std::string& str1, std::string& str2){
  bool found=false;
  unsigned short result=0;

  std::map<unsigned short, std::pair<std::string, std::string > >::iterator it; 
  for(it=_stringPairMap.begin(); it!=_stringPairMap.end(); ++it){
    std::pair<std::string, std::string >& currentPair = it->second;
    if(str1==currentPair.first && str2==currentPair.second){
      found=true;
      result=it->first;
      break;
    }   
  }

  if(!found){
    std::pair<std::string, std::string > newPair(str1, str2);
    _stringPairMapCounter++;
    result=_stringPairMapCounter;        
    _stringPairMap[result]=newPair;
  }

  return result;
}

std::pair<std::string, std::string > IdStringMapRegistry::stringPair(unsigned short id){
  if(id > _stringPairMapCounter){
    Alert << "id=" << id << " > _stringPairMapCounter=" <<  _stringPairMapCounter << endmsg;
    exit(1); 
  }
  return _stringPairMap.at(id); 
}


unsigned short IdStringMapRegistry::stringId(std::string& str1){
  bool found=false;
  unsigned short result=0;

  std::map<unsigned short, std::string >::iterator it; 
  for(it=_stringMap.begin(); it!=_stringMap.end(); ++it){
    std::string currentString = it->second;
    if(str1==currentString){
      found=true;
      result=it->first;
      break;
    }   
  }

  if(!found){
    _stringMapCounter++;
    result=_stringMapCounter;        
    _stringMap[result]=str1;
  }

  return result;
}

std::string IdStringMapRegistry::getString(unsigned short id){
  if(id > _stringMapCounter){
    Alert << "id=" << id << " > _stringMapCounter=" <<  _stringMapCounter << endmsg;
    exit(1); 
  }
  return _stringMap.at(id); 
}

