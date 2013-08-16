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

// WaveContribution class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include <memory>
#include "PwaUtils/FitParamsBase.hh"
#include "Minuit2/MnUserParameters.h"

struct calcContributionData {
   calcContributionData(std::string contribName, std::vector<std::string>& contribZeroAmpVec) :
     _contribName(contribName)
    ,_contribZeroAmpVec(contribZeroAmpVec)
  {}

  std::string _contribName;
  std::vector<std::string> _contribZeroAmpVec;

  virtual bool operator==(const calcContributionData& compare) const {
    bool result=false;
    if ( _contribName==compare._contribName) result=true;
    return result;
  }

 virtual bool operator<(const calcContributionData& compare) const {
   bool result=false;
   if(_contribName < compare._contribName) result=true; 
    return result; 
  }  
};

class AbsLh;
class fitParams;
class EvtData;
class PwaCovMatrix;
class AbsEnv;


class WaveContribution{

   public:
    WaveContribution(AbsEnv* theEnv, std::shared_ptr<AbsLh> theLh, fitParams& theFitParams);
    WaveContribution(AbsEnv* theEnv, std::shared_ptr<AbsLh> theLh, fitParams& theFitParams, 
		     std::shared_ptr<PwaCovMatrix> thePwaCovMatrix);

    std::pair<double,double> CalcContribution();
    std::vector<std::pair<std::string,std::pair<double,double>>> CalcSingleContributions();
    double CalcError(double result);

   private:
    AbsEnv* _absEnv;
    bool _calcError;
    std::shared_ptr<AbsLh> _theLh;
    std::shared_ptr<PwaCovMatrix> _thePwaCovMatrix;
    FitParamsBase _theFitParamsBase;
    ROOT::Minuit2::MnUserParameters _theMnUserParameters;
    fitParams* _theFitParamsOriginal;
    std::vector<EvtData*> _MCDataList;

    double CalcContribution(fitParams& theFitParams);
};
