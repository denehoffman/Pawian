//************************************************************************//
//                                                                        //
//  Copyright 2014 Maya Elspass (maya.elspass@ruhr-uni-bochum.de)         //
//                 Bertram Kopf (bertram@ep1.rub.de)                      //
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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include "MinFunctions/AbsFcn.hh"
#include "MinFunctions/AbsPawianMinimizer.hh"
#include "Minuit2/FCNGradientBase.h"

using namespace ROOT::Minuit2;

class AdamMinimizer : public AbsPawianMinimizer<FCNGradientBase>{

public:

    // create/copy/destroy:

    ///Constructor 
      AdamMinimizer(std::shared_ptr<AbsFcn<FCNGradientBase>> theAbsFcnPtr, std::shared_ptr<AbsPawianParameters> upar);
      
      
  ///Destructor 
      virtual ~AdamMinimizer();
  
  //Funktionen
  virtual std::string type() {return "AdamMinimizer";};
  virtual void minimize();
  virtual void printFitResult(double evtWeightSumData);
  
  // Getters:
      
protected:
  const unsigned int _max_iterations;
  std::vector<double> _s;
  std::vector<double> _v;
  double _learning_rate;
  double _initial_lr;
  double _decay_rate;
  unsigned int _iterations;
  std::shared_ptr<AbsPawianParameters> _currentPawianParams;
  double _bestLH;
  std::shared_ptr<AbsPawianParameters> _bestPawianParams;
  unsigned int _noItersWoImprovement;
  double _patience;

  void updateParameters(std::shared_ptr<AbsPawianParameters> pawianParams, std::vector<double>& gradients, std::vector<double>& s, std::vector<double>& v, int t);

 private:
};
