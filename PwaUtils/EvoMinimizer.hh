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
#include "AbsFcn.hh"
#include "Minuit2/MnUserParameters.h"

#include <boost/random/normal_distribution.hpp>


class EvoMinimizer
{
public:
   EvoMinimizer(AbsFcn& theAbsFcn, MnUserParameters upar, int iterations, int population);
   
   std::vector<double> Minimize();

private:
   int _population;
   int _iterations;
   AbsFcn* _theAbsFcn;
   MnUserParameters _bestParamsGlobal;
   MnUserParameters _tmpParams;

   void ShuffleParams();
   void AdjustSigma(double factor);

   static const char* evologo;
};
