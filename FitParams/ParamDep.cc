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

#include "FitParams/ParamDep.hh"
#include "FitParams/AbsPawianParameters.hh"

// ParamDep::ParamDep(std::vector<std::string> targetParameter, std::shared_ptr<AbsPawianParameters> params){
//   _params=params;
//   std::vector<std::string>::iterator it;
//   for(it=targetParameter.begin(); it!=targetParameter.end(); ++it){
//     _targetNames.push_back(*it);
//     _idsTarget.push_back(params->Index(*it));      
//   }  
// }

void ParamDep::Fill(std::vector<std::string> targetParameter, std::shared_ptr<AbsPawianParameters> params){
 _params=params;
  std::vector<std::string>::iterator it;
  for(it=targetParameter.begin(); it!=targetParameter.end(); ++it){
    _targetNames.push_back(*it);
    _idsTarget.push_back(params->Index(*it));      
  }
}
// void ParamDep::Fill(std::vector<std::string> targetParameter, std::istringstream& configLine, 
// 	  std::shared_ptr<AbsPawianParameters> params){
//   std::vector<std::string>::iterator it;
//   for(it=targetParameter.begin(); it!=targetParameter.end(); ++it){
//       _targetNames.push_back(*it);
//       _idsTarget.push_back(params->Index(*it));      
//   }
//   FillDerived(configLine, params);
// }
