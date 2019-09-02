//************************************************************************//
//                                                                        //
//  Copyright 2019 Bertram Kopf (bertram@ep1.rub.de)                      //
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

// ScatteringContribution class definition file. -*- C++ -*-
// Copyright 2019 Bertram Kopf

#include "qaErrorExtract/ScatteringContribution.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/EvtDataScatteringList.hh"
#include "FitParams/PwaCovMatrix.hh"
#include "ConfigParser/ParserBase.hh"
#include "ErrLogger/ErrLogger.hh"
#include "pipiScatteringUtils/pipiScatteringBaseLh.hh"
#include "Utils/IdStringMapRegistry.hh"

#include <iostream>

ScatteringContribution::ScatteringContribution(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> theFitParams,
				   std::shared_ptr<PwaCovMatrix> thePwaCovMatrix) :
  AbsCovMatErrorExtract(theLh, theFitParams, thePwaCovMatrix)
{
  _dataList=_theLh->getDataVec();
   initScatContrib();
 }

void ScatteringContribution::CalcError(){
  std::vector<EvtData*>::iterator itDataList;
   for( itDataList=_dataList.begin(); itDataList!=_dataList.end(); ++itDataList){
     CalcErrorFromData((*itDataList));
   }
}

double ScatteringContribution::CalcError(double result, std::shared_ptr<AbsPawianParameters> currentParameters){
  //dummy
   double resultErr=0;
   return resultErr;
}

void ScatteringContribution::CalcErrorFromData(EvtData* theData){

  double currentMass=theData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::M_PIPISCAT_NAME)); 

  double result=_pipiScatteringBaseLh->calcFitVal(theData, _fitParamsOriginal);

  double resultErr=0;
  double stepSize = 0.00001;
  std::map<std::string, double > derivatives;
  
  unsigned int nPar = _fitParamsOriginal->Params().size();
  std::shared_ptr<AbsPawianParameters> newFitParams = std::shared_ptr<AbsPawianParameters>(_fitParamsOriginal->Clone());
  
  for(unsigned int i=0; i<nPar; i++){
    double parOrig = _fitParamsOriginal->Value(i);
    std::string parName = _fitParamsOriginal->GetName(i);
    
    newFitParams->SetValue(i, parOrig + stepSize);
    
    double newVal = _pipiScatteringBaseLh->calcFitVal(theData, newFitParams);
    double newDerivative = (newVal - result) / stepSize;
    derivatives[parName] = newDerivative;

       newFitParams->SetValue(i, parOrig);
   }

  for(unsigned int i=0; i<nPar; i++){
    if(_theLh->CheckDoubleEquality(_fitParamsOriginal->Value(i), 0) ||
       _theLh->CheckDoubleEquality(_thePwaCovMatrix->GetElement(_fitParamsOriginal->GetName(i),
								_fitParamsOriginal->GetName(i)), 0))
      continue;

    for(unsigned int j=0; j<nPar; j++){
      if(_theLh->CheckDoubleEquality(_fitParamsOriginal->Value(j), 0))
	continue;

      std::string name1 = _fitParamsOriginal->GetName(i);
      std::string name2 = _fitParamsOriginal->GetName(j);

      resultErr += derivatives[name1] *
	_thePwaCovMatrix->GetElement(name1, name2) *
	derivatives[name2];
    }
  }
  theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FITERR_PIPISCAT_NAME))=resultErr;
  InfoMsg << "calculated error for:\t" << "m: " << currentMass << "\tval: " << result << "\terr: " << resultErr <<endmsg;
}

void ScatteringContribution::initScatContrib(){
  _pipiScatteringBaseLh=std::dynamic_pointer_cast<pipiScatteringBaseLh>(_theLh);

}
