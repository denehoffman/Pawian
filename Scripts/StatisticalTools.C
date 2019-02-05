//************************************************************************//
//                                                                        //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)                      //
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

//#include "Root.h"
#include <cmath>
#include <iostream>

class StatisticalTools
{
 public:
  static double sigmaFromLikelihoodRatio(unsigned double diffLH, unsigned int diffParams){
    return ROOT::Math::chisquared_quantile_c(TMath::Prob(2.*diffLH, diffParams), 1);
  }

  static double sigmaFromLikelihoodRatio(double LH1, int ndf1, double LH2, int ndf2){
    double diffLH=fabs(LH1-LH2);
    int diffParams=ndf1-ndf2;
    if(diffParams<0){
      cout << "ndf2 > ndf1 !!!!\n" 
           << "this is not supported!!!" << endl;
      return -10000000000000.;
    }
    return sigmaFromLikelihoodRatio(diffLH, diffParams);    
}
};
