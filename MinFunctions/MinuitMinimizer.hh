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

#pragma once

#include <vector>
#include <memory>

#include "MinFunctions/AbsFcn.hh"
#include "MinFunctions/AbsPawianMinimizer.hh"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"

#include <boost/random/normal_distribution.hpp>

using namespace ROOT::Minuit2;

class MinuitMinimizer : public AbsPawianMinimizer
{
public:
  MinuitMinimizer(std::shared_ptr<AbsFcn> theAbsFcnPtr, std::shared_ptr<AbsPawianParameters> upar);

  virtual std::string type() {return "MinuitMinimizer";};
  virtual void minimize();
  virtual void printFitResult(double evtWeightSumData);
  virtual void dumpFitResult();

protected:

  std::shared_ptr<FunctionMinimum> _mnFunctionMinimumFinalPtr;
  std::shared_ptr<MnUserParameters> _startMnUserParametersPtr;

private:
};
