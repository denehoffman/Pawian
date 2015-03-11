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
#include "FitParams/FitParColBase.hh"
#include "FitParams/PwaCovMatrix.hh"
#include "FitParams/MnPawianParameters.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ConfigParser/ParserBase.hh"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserCovariance.h"

MinuitMinimizer::MinuitMinimizer(std::shared_ptr<AbsFcn> theAbsFcnPtr, std::shared_ptr<AbsPawianParameters> upar) :
  AbsPawianMinimizer(theAbsFcnPtr, upar)
{
}


// Minimization takes place here
void MinuitMinimizer::minimize(){
  MnMigrad migrad(*_absFcn, _startPawianParams->mnUserParameters());
  Info <<"start migrad "<< endmsg;
  FunctionMinimum currentFunctionMinimum = migrad(0, GlobalEnv::instance()->parser()->tolerance());

  if(currentFunctionMinimum.IsValid()){
    //     return funcMin;
    _minimumReached=true;
    _mnFunctionMinimumFinalPtr=std::shared_ptr<FunctionMinimum>(new FunctionMinimum(currentFunctionMinimum));
    _bestPawianParams=std::shared_ptr<AbsPawianParameters>(new MnPawianParameters(_mnFunctionMinimumFinalPtr->UserParameters()));
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
       MnUserParameters newParams = _startPawianParams->mnUserParameters();
       for(unsigned int i=0; i< currentFunctionMinimum.UserParameters().Params().size();i++){
          newParams.SetValue(i, currentFunctionMinimum.UserParameters().Params().at(i));
       }
       MnMigrad migrad2(*_absFcn, newParams, MnStrategy(2));
       currentFunctionMinimum = migrad2(0, GlobalEnv::instance()->parser()->tolerance());
    }
    else{
       MnUserParameters newParams = _startPawianParams->mnUserParameters();
       for(unsigned int i=0; i< currentFunctionMinimum.UserParameters().Params().size();i++){
          newParams.SetValue(i, currentFunctionMinimum.UserParameters().Params().at(i));
          newParams.SetError(i, currentFunctionMinimum.UserParameters().Errors().at(i));
       }
       MnMigrad migrad2(*_absFcn, newParams, MnStrategy(2));
       currentFunctionMinimum = migrad2(0, GlobalEnv::instance()->parser()->tolerance());
    }

    if(currentFunctionMinimum.IsValid()){
       break;
    }
  }
  _minimumReached=true;
  _mnFunctionMinimumFinalPtr=std::shared_ptr<FunctionMinimum>(new FunctionMinimum(currentFunctionMinimum));
  _bestPawianParams=std::shared_ptr<AbsPawianParameters>(new MnPawianParameters(_mnFunctionMinimumFinalPtr->UserParameters()));
}


void MinuitMinimizer::printFitResult(double evtWeightSumData){
  if(!_minimumReached){
    Alert << "minimum has not been reached!!!" << endmsg;
    exit(1);
  }

    double theLh = _mnFunctionMinimumFinalPtr->Fval();

    Info << "\n\n********************** Final fit parameters *************************\n";
    Info << "\n" << _mnFunctionMinimumFinalPtr->UserParameters() << "\n";
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
}

void MinuitMinimizer::dumpFitResult(){
  fitParCol finalFitParams=_absFcn->defaultFitValParms();
  std::vector<double> finalParamVec=_bestPawianParams->Params();
  GlobalEnv::instance()->fitParColBase()->getFitParamVal(finalParamVec, finalFitParams);

  fitParCol finalFitErrs=_absFcn->defaultFitErrParms();
  std::vector<double> finalParamErrorVec=_bestPawianParams->Errors();  
  GlobalEnv::instance()->fitParColBase()->getFitParamVal(finalParamErrorVec, finalFitErrs);

  std::ostringstream finalResultname;

  std::string outputFileNameSuffix= GlobalEnv::instance()->outputFileNameSuffix();
  finalResultname << "finalResult" << outputFileNameSuffix << ".dat";
  
  std::ofstream theStream ( finalResultname.str().c_str() );
  GlobalEnv::instance()->fitParColBase()->dumpParams(theStream, finalFitParams, finalFitErrs);


  //dump covariance matrix
  MnUserCovariance theCovMatrix = _mnFunctionMinimumFinalPtr->UserCovariance();
  std::ostringstream serializationFileName;
  serializationFileName << "serializedOutput" << GlobalEnv::instance()->outputFileNameSuffix() << ".dat";
  std::ofstream serializationStream(serializationFileName.str().c_str());
  boost::archive::text_oarchive boostOutputArchive(serializationStream);
  
  if(_mnFunctionMinimumFinalPtr->HasValidCovariance()){
    const PwaCovMatrix thePwaCovMatrix(theCovMatrix, _mnFunctionMinimumFinalPtr->UserParameters(), finalFitParams);
    boostOutputArchive << thePwaCovMatrix;
  }
}
