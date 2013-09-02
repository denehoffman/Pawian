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
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/AbsFcn.hh"
#include "Minuit2/MnUserParameters.h"

#include <boost/random/normal_distribution.hpp>



class EvoMinimizer
{
public:
   EvoMinimizer(AbsFcn& theAbsFcn, MnUserParameters upar, int iterations,
		int population, std::string suffix="");

   std::vector<double> Minimize();

private:
   int _population;
   int _iterations;
   AbsFcn* _theAbsFcn;
   std::shared_ptr<FitParamsBase> _fitParamBase;
   fitParams _currentBestParams;
   fitParams _defaultFitErrParms;

   std::string _currentResultFileName;
   MnUserParameters _bestParamsGlobal;
   MnUserParameters _bestParamsIteration;
   MnUserParameters _tmpParams;
   MnUserParameters _iterationParamBackup;

   void ShuffleParams();
   void AdjustSigma(double factor, int numimprovements);

   static const char* evologo;
   static const double DECREASESIGMAFACTOR;
   static const double INCREASESIGMAFACTOR;
   static const double DECREASELOWTHRESH;
   static const double INCREASEHIGHTHRESH;
   static const double LHSPREADEXIT;
};
