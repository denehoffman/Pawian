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

// ParameterDependency class definition file. -*- C++ -*-
// Copyright 2014 Julian Pychy

#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <memory>

class AbsPawianParameters;

class ParamDep
{
public:
  ParamDep(){;}
  ~ParamDep(){;}

  void Fill(std::vector<std::string> targetParameter, std::shared_ptr<AbsPawianParameters> params);
  //  void Fill(std::vector<std::string> targetParameter, std::istringstream& configLine, std::shared_ptr<AbsPawianParameters> params);
  virtual void FillDerived(std::istringstream& configLine) = 0;
  virtual void Apply(std::shared_ptr<AbsPawianParameters> params) = 0;
  std::vector<unsigned int> const targetIds() const {return _idsTarget;}
  std::vector<std::string> const targetNames() {return _targetNames;}
  std::vector<unsigned int> const idRefs() {return _idRefs;}  

protected:
  std::vector<unsigned int> _idsTarget;
  std::vector<std::string> _targetNames;
  std::vector<unsigned int> _idRefs;
  std::shared_ptr<AbsPawianParameters> _params;
};

