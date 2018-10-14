//************************************************************************//
//                                                                        //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)                      //
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

// Copyright 2014 Julian Pychy, 2018 Malte Albrecht


#include "FitParams/ParamDepEqual.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "ErrLogger/ErrLogger.hh"

ParamDepEqual::ParamDepEqual(std::istringstream& configLine, std::shared_ptr<AbsPawianParameters> params){
  std::string refParameter;
  std::vector<std::string> targetParameterVec;
  configLine >> refParameter;
  _idRef = params->Index(refParameter);
  _idRefs.push_back(_idRef);
//  InfoMsg << "ParamDepEquals called with refParamId=" << _idRefs.at(0) << " and name=" << refParameter << endmsg;
  
  while(configLine) { 
    std::string parNameTarget;
    configLine >> parNameTarget;
    if(parNameTarget!="") {
      targetParameterVec.push_back(parNameTarget);
//      InfoMsg << "Target params are: " << parNameTarget << endmsg;
    }
  }

  Fill(targetParameterVec, params);
  FillDerived(configLine);     
}

ParamDepEqual::~ParamDepEqual(){
}

void ParamDepEqual::FillDerived(std::istringstream& configLine){
} 

void ParamDepEqual::Apply(std::shared_ptr<AbsPawianParameters> params){

  for(int it=0; it< _idsTarget.size(); it++) {
    params->SetValue(_idsTarget.at(it), params->Value(_idRefs.at(0)));
//    InfoMsg << "Apply: setting _idsTarget.at(" << it << ")=" << _idsTarget.at(it) << " to Value(" << _idRefs.at(0) << ")=" << params->Value(_idRefs.at(0)) << endmsg;
  }
}
