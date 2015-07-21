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

// Copyright 2014 Julian Pychy


#include "FitParams/ParamDepEqual.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "ErrLogger/ErrLogger.hh"

ParamDepEqual::ParamDepEqual(std::istringstream& configLine, std::shared_ptr<AbsPawianParameters> params){
  std::string targetParameter;
  configLine >> targetParameter;

  std::vector<std::string> targetParameterVec;
  targetParameterVec.push_back(targetParameter);
  Fill(targetParameterVec, params);
  FillDerived(configLine);     
}

ParamDepEqual::~ParamDepEqual(){
}

void ParamDepEqual::FillDerived(std::istringstream& configLine){
  std::string parNameRef;
  configLine >> parNameRef;
  _idRef = _params->Index(parNameRef);
  _idRefs.push_back(_idRef);
}

void ParamDepEqual::Apply(std::shared_ptr<AbsPawianParameters> params){
  params->SetValue(_idsTarget.at(0), params->Value(_idRef));
}





