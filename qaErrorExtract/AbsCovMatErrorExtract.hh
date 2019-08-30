//************************************************************************//
//                                                                        //
//  Copyright 2019 Bertram Kopf (bertram@ep1.rub.de)                      //
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

// AbsCovMatErrorExtract class definition file. -*- C++ -*-
// Copyright 2019 Bertram Kopf

#pragma once

#include <memory>
//#include "FitParams/FitParColBase.hh"
#include "FitParams/AbsPawianParameters.hh"

class AbsLh;
class EvtData;
class PwaCovMatrix;


class AbsCovMatErrorExtract{
  
public:
  AbsCovMatErrorExtract(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> theFitParams);
  AbsCovMatErrorExtract(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> theFitParams,
		     std::shared_ptr<PwaCovMatrix> thePwaCovMatrix);

  virtual double CalcError(double result, std::shared_ptr<AbsPawianParameters> currentParameters)=0;

protected:
  bool _calcError;
  std::shared_ptr<AbsLh> _theLh;
  std::shared_ptr<PwaCovMatrix> _thePwaCovMatrix;
  std::shared_ptr<AbsPawianParameters> _fitParamsOriginal;
  
  std::vector<EvtData*> _MCDataList;

private:
};
