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

#include <iostream>
#include <fstream>
#include <boost/random.hpp>

#include "MinFunctions/MinuitMinimizer.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "FitParams/PwaCovMatrix.hh"
#include "FitParams/MnPawianParameters.hh"
#include "FitParams/ParamDepHandler.hh"

#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ConfigParser/ParserBase.hh"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserCovariance.h"

MinuitMinimizer::MinuitMinimizer(std::shared_ptr<AbsFcn> theAbsFcnPtr, std::shared_ptr<AbsPawianParameters> upar) :
  AbsPawianMinimizer(theAbsFcnPtr, upar)
  ,_startMnUserParametersPtr(_startPawianParams->mnUserParametersPtr())
{
}


// Minimization takes place here
void MinuitMinimizer::minimize(){
  //  MnMigrad migrad(*_absFcn, _startPawianParams->mnUserParameters());
  //  MnUserParameters startMnUserP(*_startMnUserParametersPtr);
 MnMigrad migrad(*_absFcn, *_startMnUserParametersPtr);
 //    MnMigrad migrad(*_absFcn, startMnUserP);
  Info <<"start migrad "<< endmsg;
  FunctionMinimum currentFunctionMinimum = migrad(0, GlobalEnv::instance()->parser()->tolerance());

  if(currentFunctionMinimum.IsValid()){
    //     return funcMin;
    _minimumReached=true;
    _mnFunctionMinimumFinalPtr=std::shared_ptr<FunctionMinimum>(new FunctionMinimum(currentFunctionMinimum));
    //    _bestPawianParams=std::shared_ptr<AbsPawianParameters>(new MnPawianParameters(_mnFunctionMinimumFinalPtr->UserParameters()));
    _bestPawianParams->SetAllValues(_mnFunctionMinimumFinalPtr->UserParameters().Params());
    ParamDepHandler::instance()->ApplyDependencies(_bestPawianParams);
    return;
  }

  // Two more tries to get a valid result unsing strategy 2
  for(int j=0; j<2; j++){
     Warning <<"FM is invalid, try with strategy = 2."<< endmsg;

     // Check minimum covariance matrix
     bool badCovarianceDiagonal=false;
     if(currentFunctionMinimum.HasCovariance()){
        badCovarianceDiagonal = !PwaCovMatrix::DiagonalIsValid(currentFunctionMinimum.UserCovariance());
     }

     if(badCovarianceDiagonal){
       Warning << "Using default errors" << endmsg;
       std::shared_ptr<MnUserParameters> newMnUserParams = _startPawianParams->mnUserParametersPtr();
       for(unsigned int i=0; i< currentFunctionMinimum.UserParameters().Params().size();i++){
          newMnUserParams->SetValue(i, currentFunctionMinimum.UserParameters().Params().at(i));
       }
       MnMigrad migrad2(*_absFcn, *newMnUserParams, MnStrategy(2));
       currentFunctionMinimum = migrad2(0, GlobalEnv::instance()->parser()->tolerance());
    }
    else{
       std::shared_ptr<MnUserParameters> newMnUserParams = _startPawianParams->mnUserParametersPtr();
       for(unsigned int i=0; i< currentFunctionMinimum.UserParameters().Params().size();i++){
          newMnUserParams->SetValue(i, currentFunctionMinimum.UserParameters().Params().at(i));
          newMnUserParams->SetError(i, currentFunctionMinimum.UserParameters().Errors().at(i));
       }
       MnMigrad migrad2(*_absFcn, *newMnUserParams, MnStrategy(2));
       currentFunctionMinimum = migrad2(0, GlobalEnv::instance()->parser()->tolerance());
    }

    if(currentFunctionMinimum.IsValid()){
       break;
    }
  }
  _minimumReached=true;
  _mnFunctionMinimumFinalPtr=std::shared_ptr<FunctionMinimum>(new FunctionMinimum(currentFunctionMinimum));
  //  _bestPawianParams=std::shared_ptr<AbsPawianParameters>(new MnPawianParameters(_mnFunctionMinimumFinalPtr->UserParameters()));
  _bestPawianParams->SetAllValues(_mnFunctionMinimumFinalPtr->UserParameters().Params());
  ParamDepHandler::instance()->ApplyDependencies(_bestPawianParams);  
}


void MinuitMinimizer::printFitResult(double evtWeightSumData){
  if(!_minimumReached){
    Alert << "minimum has not been reached!!!" << endmsg;
    exit(1);
  }

    double theLh = _mnFunctionMinimumFinalPtr->Fval();

    Info << "\n\n********************** Final fit parameters *************************\n";
    _bestPawianParams->print(std::cout, true);
    Info << "\n\n**************** Minuit FunctionMinimum information ******************" << endmsg;
    if(_mnFunctionMinimumFinalPtr->IsValid())             Info << "\n Function minimum is valid.\n";
    else                          Info << "\n WARNING: Function minimum is invalid!" << endmsg;
    if(_mnFunctionMinimumFinalPtr->HasValidCovariance())  Info << "\n Covariance matrix is valid." << endmsg;
    else                          Info << "\n WARNING: Covariance matrix is invalid!" << endmsg;
    Info <<"\n Final LH: "<< std::setprecision(10) << theLh << "\n" << endmsg;
    Info <<" # of function calls: " << _mnFunctionMinimumFinalPtr->NFcn() << endmsg;
    Info <<" minimum edm: " << std::setprecision(10) << _mnFunctionMinimumFinalPtr->Edm()<<endmsg;
    if(!_mnFunctionMinimumFinalPtr->HasValidParameters()) Info << " hasValidParameters() returned FALSE" << endmsg;
    if(!_mnFunctionMinimumFinalPtr->HasAccurateCovar())   Info << " hasAccurateCovar() returned FALSE" << endmsg;
    if(!_mnFunctionMinimumFinalPtr->HasPosDefCovar()){    Info << " hasPosDefCovar() returned FALSE" << endmsg;
                                  if(_mnFunctionMinimumFinalPtr->HasMadePosDefCovar()) Info << " hasMadePosDefCovar() returned TRUE" << endmsg;
    }
    if(!_mnFunctionMinimumFinalPtr->HasCovariance())      Info << " hasCovariance() returned FALSE" << endmsg;
    if(_mnFunctionMinimumFinalPtr->HasReachedCallLimit()) Info << " hasReachedCallLimit() returned TRUE" << endmsg;
    if(_mnFunctionMinimumFinalPtr->IsAboveMaxEdm())       Info << " isAboveMaxEdm() returned TRUE" << endmsg;
    if(_mnFunctionMinimumFinalPtr->HesseFailed())         Info << " hesseFailed() returned TRUE\n" << endmsg;


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // calculate AIC, BIC criteria and output selected wave contrib
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    unsigned int noOfFreeFitParams=_bestPawianParams->VariableParameters();

    double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
    double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
    double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/(evtWeightSumData-noOfFreeFitParams-1);
    Info << "noOfFreeFitParams:\t" <<noOfFreeFitParams;
    Info << "evtWeightSumData:\t" <<evtWeightSumData;
    Info << "BIC:\t" << BICcriterion << endmsg;
    Info << "AIC:\t" << AICcriterion << endmsg;
    Info << "AICc:\t" << AICccriterion << endmsg;
}

void MinuitMinimizer::dumpFitResult(){

  std::ostringstream finalResultname;
  std::string outputFileNameSuffix= GlobalEnv::instance()->outputFileNameSuffix();
  finalResultname << "finalResult" << outputFileNameSuffix << ".dat";
  
  std::ofstream theStream ( finalResultname.str().c_str() );
  _bestPawianParams->print(theStream);

  //dump covariance matrix
  MnUserCovariance theCovMatrix = _mnFunctionMinimumFinalPtr->UserCovariance();
  std::ostringstream serializationFileName;
  serializationFileName << "serializedOutput" << GlobalEnv::instance()->outputFileNameSuffix() << ".dat";
  std::ofstream serializationStream(serializationFileName.str().c_str());
  boost::archive::text_oarchive boostOutputArchive(serializationStream);
  
  if(_mnFunctionMinimumFinalPtr->HasValidCovariance()){
    const PwaCovMatrix thePwaCovMatrix(theCovMatrix, _mnFunctionMinimumFinalPtr->UserParameters());
    boostOutputArchive << thePwaCovMatrix;
  }
}
