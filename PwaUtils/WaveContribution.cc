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


#include "PwaUtils/WaveContribution.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/PwaCovMatrix.hh"
#include "ErrLogger/ErrLogger.hh"


#include <iostream>

WaveContribution::WaveContribution(std::shared_ptr<AbsLh> theLh, fitParams& theFitParams) :
      _calcError(false)
    , _theLh(theLh)
    , _theFitParamsOriginal(&theFitParams)
{
   _MCDataList=_theLh->getMcVec();
   _theFitParamsBase.setMnUsrParams(_theMnUserParameters, *_theFitParamsOriginal, *_theFitParamsOriginal);
}



WaveContribution::WaveContribution(std::shared_ptr<AbsLh> theLh, fitParams& theFitParams, 
				   std::shared_ptr<PwaCovMatrix> thePwaCovMatrix) :
    _calcError(true)
   ,_theLh(theLh)
   ,_thePwaCovMatrix(thePwaCovMatrix)
   , _theFitParamsOriginal(&theFitParams)
{
   _MCDataList=_theLh->getMcVec();
   _theFitParamsBase.setMnUsrParams(_theMnUserParameters, *_theFitParamsOriginal, *_theFitParamsOriginal);
}





double WaveContribution::CalcContribution(fitParams& theFitParams){

   double result=0;

   for(auto it=_MCDataList.begin(); it!=_MCDataList.end(); ++it){
      result +=_theLh->calcEvtIntensity( (*it), theFitParams);   
   } 

   return result;
}



std::pair<double,double> WaveContribution::CalcContribution(){

   double result=CalcContribution(*_theFitParamsOriginal);
   double resultErr=0;

   if(!_calcError)
      return std::pair<double,double>(result, resultErr);

   double stepSize = 0.0001;
   std::map< std::string, double > derivatives;

   unsigned int nPar = _theMnUserParameters.Params().size();
   
   for(unsigned int i=0; i<nPar; i++){
      double parOrig = _theMnUserParameters.Value(i);
      std::string parName = _theMnUserParameters.GetName(i);

      _theMnUserParameters.SetValue(i, parOrig + stepSize);

      fitParams newFitParams = *_theFitParamsOriginal;
      _theFitParamsBase.getFitParamVal(_theMnUserParameters.Params(), newFitParams);
      _theLh->updateFitParams(newFitParams);

      double newContribution = CalcContribution(newFitParams);
      double newDerivative = (newContribution - result) / stepSize;
      derivatives[parName] = newDerivative;

      _theMnUserParameters.SetValue(i, parOrig);
   }

   for(unsigned int i=0; i<nPar; i++){

      if(_theLh->CheckDoubleEquality(_theMnUserParameters.Value(i), 0) ||
	 _theLh->CheckDoubleEquality(_thePwaCovMatrix->GetElement(_theMnUserParameters.GetName(i),
       								  _theMnUserParameters.GetName(i)), 0))
	 continue;
      
      Info << "Param used in contribution error calculation: " << _theMnUserParameters.GetName(i) << endmsg;

      for(unsigned int j=0; j<nPar; j++){

	 if(_theLh->CheckDoubleEquality(_theMnUserParameters.Value(j), 0))
	    continue;

         std::string name1 = _theMnUserParameters.GetName(i);
         std::string name2 = _theMnUserParameters.GetName(j);

	 resultErr += derivatives[name1] * 
	    _thePwaCovMatrix->GetElement(name1, name2) * 
	    derivatives[name2];
      }
   }
   
   resultErr = sqrt(resultErr);

   return std::pair<double,double>(result, resultErr);
}

