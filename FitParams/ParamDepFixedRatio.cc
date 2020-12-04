//************************************************************************//
//                                                                        //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)                      //
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

// Copyright 2020 Bertram Kopf


#include "FitParams/ParamDepFixedRatio.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "ErrLogger/ErrLogger.hh"

ParamDepFixedRatio::ParamDepFixedRatio(std::istringstream& configLine, std::shared_ptr<AbsPawianParameters> params) :
  _scalingFactor(0.)
{
  std::string refParameter;
  std::vector<std::string> targetParameterVec;
  configLine >> refParameter;
  _idRef = params->Index(refParameter);
  _idRefs.push_back(_idRef);
   
  std::string parNameTarget;
  configLine >> parNameTarget;
  if(parNameTarget!="") {
    targetParameterVec.push_back(parNameTarget);
    _idRef = params->Index(refParameter);
  _idRefs.push_back(_idRef);
  }

  std::string scaleStr;
  configLine >> scaleStr;
  _scalingFactor = atof(scaleStr.c_str());

  InfoMsg << "refParam: " << refParameter << " targetParam:" <<  parNameTarget << " scalingFactor: " << _scalingFactor << endmsg;
  
  Fill(targetParameterVec, params);
  FillDerived(configLine);     
}

ParamDepFixedRatio::~ParamDepFixedRatio(){
}

void ParamDepFixedRatio::FillDerived(std::istringstream& configLine){
} 

void ParamDepFixedRatio::Apply(std::shared_ptr<AbsPawianParameters> params){
  //for(size_t it=0; it< _idsTarget.size(); it++) {
  //  params->SetValue(_idsTarget.at(it), _scalingFactor*params->Value(_idRefs.at(0)));
  //}
    params->SetValue(_idsTarget.at(0), _scalingFactor*params->Value(_idRefs.at(0)));
}
