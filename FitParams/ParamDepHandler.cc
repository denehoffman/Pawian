//************************************************************************//
//                                                                        //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 Julian Pychy (julian@ep1.rub.de)                       //
//                 - Ruhr-Universität Bochum                              //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

#include "FitParams/ParamDepHandler.hh"
#include "FitParams/ParamDepEqual.hh"
#include "FitParams/ParamDepGFactorToFixFullWidth.hh"
#include "FitParams/ParamDepGFacToFixGFacsWidthMass.hh"

#include "ErrLogger/ErrLogger.hh"

ParamDepHandler* ParamDepHandler::_instance = 0;

ParamDepHandler* ParamDepHandler::instance(){
  if(_instance == 0){
     _instance = new ParamDepHandler();
  }
  return _instance;
}



void ParamDepHandler::Fill(const std::vector<std::string>& configLines, 
				      std::shared_ptr<AbsPawianParameters> params){

  for(auto it = configLines.begin(); it!=configLines.end(); ++it){
     std::istringstream currentStream(*it);
     
     std::shared_ptr<ParamDep> newDependency;
     std::string targetParameter;
     std::string type;
     currentStream >> type;
     
     if(type == "equals"){
       newDependency = std::shared_ptr<ParamDep>(new ParamDepEqual(currentStream, params));
     }
     else if(type == "gFactorToFixFullWidth"){
	newDependency = std::shared_ptr<ParamDep>(new ParamDepGFactorToFixFullWidth(currentStream, params));
     }
     else if(type == "gFacToFixGFacsWidthMass"){
       newDependency = std::shared_ptr<ParamDep>(new ParamDepGFacToFixGFacsWidthMass(currentStream, params));
     }
     else{
	Alert << "Dependency type not found: " << type << endmsg;
	exit(0);
     }
     
     Info << "Adding parameter dependency for " << targetParameter << " type "
	  << type << " arguments " << currentStream.str() << endmsg;

     std::vector<std::string> targetParameterVec=newDependency->targetNames();
     std::vector<std::string>::iterator itParam;
     for(itParam= targetParameterVec.begin(); itParam!= targetParameterVec.end(); ++itParam){
       std::vector<std::string>::iterator foundParam;
       foundParam=std::find(_dependentParameterNames.begin(), _dependentParameterNames.end(), (*itParam));
       if(foundParam == _dependentParameterNames.end()) _dependentParameterNames.push_back(*itParam);
     }
     _dependencies.push_back(newDependency);
  }
}

void ParamDepHandler::ApplyDependencies(std::shared_ptr<AbsPawianParameters> params){
  for(auto it = _dependencies.begin(); it!= _dependencies.end(); ++it){
     (*it)->Apply(params);
  }
}

std::vector<unsigned int> ParamDepHandler::RefIds(std::string parName){
  std::vector<unsigned int> result;
  for(auto it = _dependencies.begin(); it!= _dependencies.end(); ++it){
    std::vector<std::string> currentNames=(*it)->targetNames();
    std::vector<std::string>::const_iterator itNames;
    for(itNames=currentNames.begin(); itNames!=currentNames.end(); ++itNames){
      if( (*itNames)==parName) result=(*it)->idRefs();
      continue;
    }
  }
  return result;
} 
