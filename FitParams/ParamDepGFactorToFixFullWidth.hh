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

// ParamDepGFactorToFixFullWidth class definition file. -*- C++ -*-
// Copyright 2014 Julian Pychy

#include "FitParams/ParamDep.hh"

#include <vector>
#include <string>
#include <memory>
#include <sstream>

class AbsPawianParameters; 

struct ChannelData
{
  int massParamId;
  int gParamId;
  double m1;
  double m2;
  
};


class ParamDepGFactorToFixFullWidth : public ParamDep
{
public:
  virtual void FillDerived(std::istringstream& configLine,
		    std::shared_ptr<AbsPawianParameters> params);
  virtual void Apply(std::shared_ptr<AbsPawianParameters> params);

private:
  ChannelData _targetChannelData;
  std::vector<ChannelData> _refChannelData;
  double _targetFullWidth;
};



