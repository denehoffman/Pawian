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

// ScatteringContribution class definition file. -*- C++ -*-
// Copyright 2019 Bertram Kopf

#pragma once

#include <memory>
#include "qaErrorExtract/AbsCovMatErrorExtract.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "pipiScatteringUtils/pipiScatteringBaseLh.hh"

class AbsLh;
class EvtData;
class PwaCovMatrix;
class pipiScatteringBaseLh;

class ScatteringContribution : public AbsCovMatErrorExtract{
  
  public:
  ScatteringContribution(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> theFitParams,
		     std::shared_ptr<PwaCovMatrix> thePwaCovMatrix);

  void CalcError();
  virtual double CalcError(double result, std::shared_ptr<AbsPawianParameters> currentParameters);

protected:
  std::vector<EvtData*> _dataList;
  std::shared_ptr<pipiScatteringBaseLh> _pipiScatteringBaseLh;

  void CalcErrorFromData(EvtData* theData);
private:
  void initScatContrib();
};
