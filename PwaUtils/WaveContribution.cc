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
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/PwaCovMatrix.hh"
#include "ConfigParser/ParserBase.hh"
#include "epemUtils/epemHist.hh"
#include "pbarpUtils/pbarpHist.hh"
#include "ErrLogger/ErrLogger.hh"


#include <iostream>

WaveContribution::WaveContribution(std::shared_ptr<AbsLh> theLh, fitParams& theFitParams) :
      _calcError(false)
    , _theLh(theLh)
    , _theFitParamsOriginal(&theFitParams)
{
   _MCDataList=_theLh->getMcVec();
}



WaveContribution::WaveContribution(std::shared_ptr<AbsLh> theLh, fitParams& theFitParams,
				   std::shared_ptr<PwaCovMatrix> thePwaCovMatrix) :
     _calcError(true)
   , _theLh(theLh)
   , _thePwaCovMatrix(thePwaCovMatrix)
   , _theFitParamsOriginal(&theFitParams)
{
   _MCDataList=_theLh->getMcVec();
}





double WaveContribution::CalcContribution(fitParams& theFitParams){

   double result=0;

   _theLh->updateFitParams(theFitParams);

   for(auto it=_MCDataList.begin(); it!=_MCDataList.end(); ++it){
      result +=_theLh->calcEvtIntensity( (*it), theFitParams);
   }

   return result;
}



std::pair<double,double> WaveContribution::CalcContribution(){

   double result=CalcContribution(*_theFitParamsOriginal);
   double resultErr=0;

   if(!_calcError){
      return std::pair<double,double>(result, resultErr);
   }
   else{
      MnUserParameters theMnUserParameters;
      _theFitParamsBase.setMnUsrParams(theMnUserParameters, *_theFitParamsOriginal, *_theFitParamsOriginal);      
      return std::pair<double,double>(result, CalcError(result, theMnUserParameters));
   }
}



std::vector<std::pair<std::string,std::pair<double,double>>> WaveContribution::CalcSingleContributions(){

   std::vector<std::pair<std::string,std::pair<double,double>>> retValues;
   std::vector<std::shared_ptr<calcContributionData> > calcContributionDataVec = GlobalEnv::instance()->Channel()->calcContributionDataVec();

   for(unsigned int i=0; i<NoOfContributions(); i++){

     std::string tmpContribName = GetContributionName(i);
     MnUserParameters currentMnUserParameters = GetParametersForContribution(i); 

     fitParams newFitParams = *_theFitParamsOriginal;
     _theFitParamsBase.getFitParamVal(currentMnUserParameters.Params(), newFitParams);
     double newContribution = CalcContribution(newFitParams); // calls updateFitParams

     if(!_calcError){
       retValues.push_back(std::pair<std::string,std::pair<double,double>>(tmpContribName, std::pair<double,double>(newContribution, 0)));
     }
     else{
       double error = CalcError(newContribution, currentMnUserParameters);
       retValues.push_back(std::pair<std::string,std::pair<double,double>>(tmpContribName, std::pair<double,double>(newContribution, error)));
     }
     if(GlobalEnv::instance()->parser()->saveContributionHistos()){
       auto hist = GlobalEnv::instance()->Channel()->CreateHistInstance(tmpContribName);
       hist->fillFromLhData(_theLh, newFitParams);
     }

     Info << "calculated contribution for " << tmpContribName << " = " << newContribution << " +- " << retValues.back().second.second << endmsg;
   }

   return retValues;
}



double WaveContribution::CalcError(double result, ROOT::Minuit2::MnUserParameters currentMnUserParameters){
   double resultErr=0;
   double stepSize = 0.0001;
   std::map< std::string, double > derivatives;

   unsigned int nPar = currentMnUserParameters.Params().size();

   for(unsigned int i=0; i<nPar; i++){
      double parOrig = currentMnUserParameters.Value(i);
      std::string parName = currentMnUserParameters.GetName(i);

      currentMnUserParameters.SetValue(i, parOrig + stepSize);

      fitParams newFitParams = *_theFitParamsOriginal;
      _theFitParamsBase.getFitParamVal(currentMnUserParameters.Params(), newFitParams);

      double newContribution = CalcContribution(newFitParams);
      double newDerivative = (newContribution - result) / stepSize;
      derivatives[parName] = newDerivative;

      currentMnUserParameters.SetValue(i, parOrig);
   }

   for(unsigned int i=0; i<nPar; i++){

      if(_theLh->CheckDoubleEquality(currentMnUserParameters.Value(i), 0) ||
	 _theLh->CheckDoubleEquality(_thePwaCovMatrix->GetElement(currentMnUserParameters.GetName(i),
       								  currentMnUserParameters.GetName(i)), 0))
	 continue;

      Info << "Param used in contribution error calculation: " << currentMnUserParameters.GetName(i) << endmsg;

      for(unsigned int j=0; j<nPar; j++){

	 if(_theLh->CheckDoubleEquality(currentMnUserParameters.Value(j), 0))
	    continue;

         std::string name1 = currentMnUserParameters.GetName(i);
         std::string name2 = currentMnUserParameters.GetName(j);

	 resultErr += derivatives[name1] *
	    _thePwaCovMatrix->GetElement(name1, name2) *
	    derivatives[name2];
      }
   }

   resultErr = sqrt(resultErr);
   return resultErr;
}



unsigned int WaveContribution::NoOfContributions(){
  return GlobalEnv::instance()->Channel()->calcContributionDataVec().size();
}



std::string WaveContribution::GetContributionName(unsigned int index){
  
  if(index >= NoOfContributions()){
     Alert << "index > NoOfContributions()" << endmsg;
  }

  return GlobalEnv::instance()->Channel()->calcContributionDataVec().at(index)->_contribName;
}



ROOT::Minuit2::MnUserParameters WaveContribution::GetParametersForContribution(unsigned int index){

  if(index >= NoOfContributions()){
     Alert << "index > NoOfContributions()" << endmsg;
  }

  MnUserParameters newMnUserParameters;
  _theFitParamsBase.setMnUsrParams(newMnUserParameters, *_theFitParamsOriginal, *_theFitParamsOriginal);
  unsigned int nPar = newMnUserParameters.Params().size();

  std::vector<std::shared_ptr<calcContributionData> > calcContributionDataVec = GlobalEnv::instance()->Channel()->calcContributionDataVec();
  std::vector<std::string> tmpZeroAmp = calcContributionDataVec.at(index)->_contribZeroAmpVec;

  for(auto itZeroAmpVec = tmpZeroAmp.begin(); itZeroAmpVec!=tmpZeroAmp.end(); ++itZeroAmpVec) {      // loop over to be zeroed amplitudes in ONE "calcContribution"-line
     for(unsigned int i=0; i<nPar; i++){  // loop over all existing fitParameters
	std::string parName = newMnUserParameters.GetName(i);
	if(parName.find(*itZeroAmpVec) != std::string::npos){
	   Info << "setting parameter to 0.0: " << parName << endmsg;
	   newMnUserParameters.SetValue(i, 0.);
	}
     }
  }

  return newMnUserParameters;
}
