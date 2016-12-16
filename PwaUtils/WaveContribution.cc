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
#include "FitParams/PwaCovMatrix.hh"
#include "ConfigParser/ParserBase.hh"
#include "epemUtils/epemHist.hh"
#include "pbarpUtils/pbarpHist.hh"
#include "ErrLogger/ErrLogger.hh"
//#include "FitParams/ParamFactory.hh"

#include <iostream>

WaveContribution::WaveContribution(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> theFitParams) :
     _calcError(false)
    , _theLh(theLh)
    , _fitParamsOriginal(theFitParams)
{
   _MCDataList=_theLh->getMcVec();
}

WaveContribution::WaveContribution(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> theFitParams,
				   std::shared_ptr<PwaCovMatrix> thePwaCovMatrix) :
     _calcError(true)
   , _theLh(theLh)
   , _thePwaCovMatrix(thePwaCovMatrix)
   , _fitParamsOriginal(theFitParams)
{
   _MCDataList=_theLh->getMcVec();
}

double WaveContribution::CalcContribution(std::shared_ptr<AbsPawianParameters> theFitParams){

   double result=0;

   _theLh->updateFitParams(theFitParams);

   for(auto it=_MCDataList.begin(); it!=_MCDataList.end(); ++it){
      result +=_theLh->calcEvtIntensity( (*it), theFitParams);
   }

   return result;
}


std::pair<double,double> WaveContribution::CalcContribution(){

   double result=CalcContribution(_fitParamsOriginal);
   double resultErr=0.;

   if(!_calcError){
      return std::pair<double,double>(result, resultErr);
   }
   else{
     std::shared_ptr<AbsPawianParameters> theParameters=std::shared_ptr<AbsPawianParameters>(_fitParamsOriginal->Clone());      
      return std::pair<double,double>(result, CalcError(result, theParameters));
   }
}



std::vector<std::pair<std::string,std::pair<double,double>>> WaveContribution::CalcSingleContributions(){

   std::vector<std::pair<std::string,std::pair<double,double>>> retValues;

   for(unsigned int i=0; i<NoOfContributions(); i++){

     std::string tmpContribName = GetContributionName(i);
     std::shared_ptr<AbsPawianParameters> currentParameters = GetParametersForContribution(i); 

     double newContribution = CalcContribution(currentParameters); // calls updateFitParams

     if(!_calcError){
       retValues.push_back(std::pair<std::string,std::pair<double,double>>(tmpContribName, std::pair<double,double>(newContribution, 0)));
     }
     else{
       double error = CalcError(newContribution, currentParameters);
       retValues.push_back(std::pair<std::string,std::pair<double,double>>(tmpContribName, std::pair<double,double>(newContribution, error)));
     }
     if(GlobalEnv::instance()->parser()->saveContributionHistos()){
       auto hist = GlobalEnv::instance()->Channel()->CreateHistInstance(tmpContribName);
       //       hist->fillFromLhData(_theLh, newFitParams);
       hist->fillFromLhData(_theLh, currentParameters);
     }

     InfoMsg << "calculated contribution for " << tmpContribName << " = " << newContribution << " +- " << retValues.back().second.second << endmsg;
   }

   return retValues;
}



double WaveContribution::CalcError(double result, std::shared_ptr<AbsPawianParameters> currentParameters){
   double resultErr=0;
   double stepSize = 0.0001;
   std::map< std::string, double > derivatives;

   unsigned int nPar = currentParameters->Params().size();
   std::shared_ptr<AbsPawianParameters> newFitParams = std::shared_ptr<AbsPawianParameters>(currentParameters->Clone());

   for(unsigned int i=0; i<nPar; i++){
      double parOrig = currentParameters->Value(i);
      std::string parName = currentParameters->GetName(i);

      newFitParams->SetValue(i, parOrig + stepSize);
 
      double newContribution = CalcContribution(newFitParams);
      double newDerivative = (newContribution - result) / stepSize;
      derivatives[parName] = newDerivative;

      newFitParams->SetValue(i, parOrig);
   }

   for(unsigned int i=0; i<nPar; i++){

      if(_theLh->CheckDoubleEquality(currentParameters->Value(i), 0) ||
	 _theLh->CheckDoubleEquality(_thePwaCovMatrix->GetElement(currentParameters->GetName(i),
       								  currentParameters->GetName(i)), 0))
	 continue;

      InfoMsg << "Param used in contribution error calculation: " << currentParameters->GetName(i) << endmsg;

      for(unsigned int j=0; j<nPar; j++){

	 if(_theLh->CheckDoubleEquality(currentParameters->Value(j), 0))
	    continue;

         std::string name1 = currentParameters->GetName(i);
         std::string name2 = currentParameters->GetName(j);

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


std::shared_ptr<AbsPawianParameters> WaveContribution::GetParametersForContribution(unsigned int index){

  if(index >= NoOfContributions()){
     Alert << "index > NoOfContributions()" << endmsg;
  }

  std::shared_ptr<AbsPawianParameters> newParameters = std::shared_ptr<AbsPawianParameters>(_fitParamsOriginal->Clone());
  unsigned int nPar = newParameters->Params().size();

  std::vector<std::shared_ptr<calcContributionData> > calcContributionDataVec = GlobalEnv::instance()->Channel()->calcContributionDataVec();
  std::vector<std::string> tmpZeroAmp = calcContributionDataVec.at(index)->_contribZeroAmpVec;

  for(auto itZeroAmpVec = tmpZeroAmp.begin(); itZeroAmpVec!=tmpZeroAmp.end(); ++itZeroAmpVec) {      // loop over to be zeroed amplitudes in ONE "calcContribution"-line
     for(unsigned int i=0; i<nPar; i++){  // loop over all existing fitParameters
	std::string parName = newParameters->GetName(i);
	if(parName.find(*itZeroAmpVec) != std::string::npos){
	   InfoMsg << "setting parameter to 0.0: " << parName << endmsg;
	   newParameters->SetValue(i, 0.);
	}
     }
  }

  return newParameters;
}
