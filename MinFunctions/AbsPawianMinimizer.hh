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

// AbsPawianMinimizer class definition file. -*- C++ -*-
// Copyright 2015 Bertram Kopf

#pragma once

//#include <iostream>
//#include <fstream>
#include <vector>
#include <memory>

#include "MinFunctions/AbsFcn.hh"
#include "FitParams/AbsPawianParameters.hh"

class AbsPawianMinimizer {

public:

  AbsPawianMinimizer(std::shared_ptr<AbsFcn> theAbsFcnPtr, std::shared_ptr<AbsPawianParameters> upar);

  ~AbsPawianMinimizer() {}

  virtual std::string type()=0;
  virtual void minimize()=0;
  virtual void printFitResult(double evtWeightSumData)=0;
  virtual void printFitResultQA(double evtWeightSumData);
  virtual void dumpFitResult();

protected:
  std::shared_ptr<AbsFcn> _absFcn;
  std::shared_ptr<AbsPawianParameters> _startPawianParams;
  std::shared_ptr<AbsPawianParameters> _bestPawianParams;
  bool _minimumReached;
  double _finalLh;  

private:

};

