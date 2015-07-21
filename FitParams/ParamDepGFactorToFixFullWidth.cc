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

#include <iostream>

#include "FitParams/ParamDepGFactorToFixFullWidth.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "ErrLogger/ErrLogger.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

ParamDepGFactorToFixFullWidth::ParamDepGFactorToFixFullWidth(std::istringstream& configLine, std::shared_ptr<AbsPawianParameters> params){
  std::string targetParameter;
  configLine >> targetParameter;

  std::vector<std::string> targetParameterVec;
  targetParameterVec.push_back(targetParameter);
  Fill(targetParameterVec, params);
  FillDerived(configLine);
}

ParamDepGFactorToFixFullWidth::~ParamDepGFactorToFixFullWidth(){
}

void ParamDepGFactorToFixFullWidth::FillDerived(std::istringstream& configLine){

  configLine >> _targetFullWidth;

  std::string massParamName;
  configLine >> massParamName;
  _targetChannelData.massParamId = _params->Index(massParamName);
  _targetChannelData.gParamId = _idsTarget.at(0);
  configLine >> _targetChannelData.m1 >> _targetChannelData.m2;

  int numRefChannels;
  configLine >> numRefChannels;

  for(int i=1; i<=numRefChannels; i++){
     ChannelData newChannelData;
     configLine >> massParamName;

     newChannelData.massParamId = _params->Index(massParamName);
     _idRefs.push_back(_params->Index(massParamName));
     std::string gParamName;
     configLine >> gParamName;

     newChannelData.gParamId = _params->Index(gParamName);
     _idRefs.push_back(_params->Index(gParamName));
     configLine >> newChannelData.m1 >> newChannelData.m2;
     _refChannelData.push_back(newChannelData);

  }
}

void ParamDepGFactorToFixFullWidth::Apply(std::shared_ptr<AbsPawianParameters> params){
  double result = 0;

  double partialWidthSum = 0;
  for(auto it = _refChannelData.begin(); it!= _refChannelData.end(); ++it){
     // Gamma_i = g*g*rho(m,m1,m2)/m
    double m = params->Value((*it).massParamId);
    partialWidthSum += (params->Value((*it).gParamId) * params->Value((*it).gParamId) * phaseSpaceFac(m, (*it).m1, (*it).m2) / m).real();
  }

  double targetPartialWidth = _targetFullWidth - partialWidthSum;
  double mTarget = params->Value(_targetChannelData.massParamId);
  
  result = (mTarget / phaseSpaceFac(mTarget, _targetChannelData.m1, _targetChannelData.m2)).real() * targetPartialWidth;

  if((targetPartialWidth < 0) || (result < 0)){
     Alert << "Invalid g-Factor value! "
	   << "Partial width sum = " <<  partialWidthSum 
	   << " dependent width = "   << targetPartialWidth 
           << " Mass = " << mTarget
	   << " g^2 = " << result << endmsg;

     for(auto it = _refChannelData.begin(); it!= _refChannelData.end(); ++it){
	Alert << "Independent Channel: g = " << params->Value((*it).gParamId) << " m = " << params->Value((*it).massParamId) << endmsg;
     }

     exit(0);
  }
  params->SetValue(_targetChannelData.gParamId, sqrt(result));
}




