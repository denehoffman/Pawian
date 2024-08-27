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
#include "Minuit2/FCNBase.h"
#include "Minuit2/FCNGradientBase.h"

#include "MinFunctions/AbsPawianMinimizer.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/GlobalEnv.hh"

template<typename T>
AbsPawianMinimizer<T>::AbsPawianMinimizer(std::shared_ptr<AbsFcn<T>> theAbsFcnPtr, std::shared_ptr<AbsPawianParameters> upar) :
  _absFcn(theAbsFcnPtr)
  ,_startPawianParams(upar)
  ,_bestPawianParams(upar->Clone())
  ,_minimumReached(false)
  ,_finalLh(10000000.)
{
}

template<typename T>
AbsPawianMinimizer<T>::~AbsPawianMinimizer(){
}


template<typename T>
void AbsPawianMinimizer<T>::dumpFitResult(){

  std::ostringstream finalResultname;

  std::string outputFileNameSuffix= GlobalEnv::instance()->outputFileNameSuffix();
  finalResultname << "finalResult" << outputFileNameSuffix << ".dat";
  
  std::ofstream theStream ( finalResultname.str().c_str() );
  _bestPawianParams->print(theStream);
}

template<typename T>
void AbsPawianMinimizer<T>::printFitResultQA(double evtWeightSumData){
  //double theLh = _mnFunctionMinimumFinalPtr->Fval();
    InfoMsg << "\n\n********************** fit parameters *************************" << endmsg;
    _bestPawianParams->print(std::cout, true);
    const std::vector<double> theParams=_bestPawianParams->Params();
    double theLH=(*_absFcn)(theParams);
    unsigned int noOfFreeFitParams=_bestPawianParams->VariableParameters();

    double BICcriterion=2.*theLH+noOfFreeFitParams*log(evtWeightSumData);
    double AICcriterion=2.*theLH+2.*noOfFreeFitParams;
    double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)
          / (evtWeightSumData-noOfFreeFitParams-1);
    InfoMsg << "NLL:\t" << theLH << endmsg;
    InfoMsg << "noOfFreeFitParams:\t" <<noOfFreeFitParams << endmsg;
    InfoMsg << "evtWeightSumData:\t" <<evtWeightSumData << endmsg;
    InfoMsg << "BIC:\t" << BICcriterion << endmsg;
    InfoMsg << "AIC:\t" << AICcriterion << endmsg;
    InfoMsg << "AICc:\t" << AICccriterion << endmsg;
}

template AbsPawianMinimizer<FCNBase>::AbsPawianMinimizer(std::shared_ptr<AbsFcn<FCNBase>>, std::shared_ptr<AbsPawianParameters>);
template AbsPawianMinimizer<FCNBase>::~AbsPawianMinimizer();
template AbsPawianMinimizer<FCNGradientBase>::AbsPawianMinimizer(std::shared_ptr<AbsFcn<FCNGradientBase>>, std::shared_ptr<AbsPawianParameters>);
template AbsPawianMinimizer<FCNGradientBase>::~AbsPawianMinimizer();
