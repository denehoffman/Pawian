//************************************************************************//
//                                                                        //
//  Copyright 2014 Maya Elspass (maya.elspass@ruhr-uni-bochum.de)         //
//                 Bertram Kopf (bertram@ep1.rub.de)                      //
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

#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <iomanip>


#include "MinFunctions/AdamMinimizer.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ConfigParser/ParserBase.hh"

AdamMinimizer::AdamMinimizer(std::shared_ptr<AbsFcn<FCNGradientBase>> theAbsFcnPtr, std::shared_ptr<AbsPawianParameters> upar) :
  AbsPawianMinimizer<FCNGradientBase>(theAbsFcnPtr, upar)
  ,_max_iterations(1000)
  ,_s(upar->Params().size(),0.)
  ,_v(upar->Params().size(),0.)
  ,_learning_rate(0.2)
  ,_iterations(0)
  ,_currentPawianParams(std::shared_ptr<AbsPawianParameters>(upar->Clone()))
  ,_bestLH(1.e20)
  ,_bestPawianParams(std::shared_ptr<AbsPawianParameters>(upar->Clone()))
  ,_noItersWoImprovement(0)
{
}

AdamMinimizer::~AdamMinimizer()
{
}

void AdamMinimizer::minimize(){
  if(_iterations==0){
    double currentLH=_absFcn->operator()(_bestPawianParams->Params());
    if(currentLH<_bestLH) _bestLH=currentLH;
    InfoMsg << "best LH: " << _bestLH << "\tnoIters: " << _iterations << "\tnoItersWoImprovement: " << _noItersWoImprovement << endmsg;    
  }
    while(_iterations <= _max_iterations && _noItersWoImprovement<100){
      if(_noItersWoImprovement>0 && _noItersWoImprovement%10 == 0){ //continue with parameters of the best fit
	_currentPawianParams=std::shared_ptr<AbsPawianParameters>(_bestPawianParams->Clone());
      }
      std::vector<double> derivatives = _absFcn->Gradient(_currentPawianParams->Params());

      updateParameters(_currentPawianParams, derivatives, _s, _v, _iterations);
      double currentLH=_absFcn->operator()(_currentPawianParams->Params());
      if(currentLH<_bestLH){
	_bestLH=currentLH;
	_bestPawianParams=std::shared_ptr<AbsPawianParameters>(_currentPawianParams->Clone());
	_noItersWoImprovement=0;
      }
      else ++_noItersWoImprovement;
      ++_iterations;
      _learning_rate *= 0.999;
      InfoMsg << "best LH: " << _bestLH << "\tnoIters: " << _iterations << "\tnoItersWoImprovement: " << _noItersWoImprovement << endmsg;
    }
}

void AdamMinimizer::updateParameters(std::shared_ptr<AbsPawianParameters> pawianParams, std::vector<double>& gradients, std::vector<double>& s, std::vector<double>& v, int t){
  double beta1=0.9;
  double beta2=0.99; 
  double epsilon=1.*pow(10., -8);

  for(unsigned int i = 0; i < pawianParams->Params().size(); ++i){
    if (pawianParams->IsFixed(i)) continue;
    _s.at(i) = beta1 * _s.at(i) + (1.0 - beta1) * gradients.at(i);
    _v.at(i) = beta2 * _v.at(i) + (1.0 - beta2) * gradients.at(i) * gradients.at(i);

    double s_hat = _s.at(i) / (1.0 - pow(beta1 , (t + 1)));
    double v_hat = _v.at(i) / (1.0 - pow(beta2 , (t + 1)));

    double newVal = pawianParams->Value(i) - _learning_rate * s_hat / (std::sqrt(v_hat) + epsilon);
    if(pawianParams->HasLimits(i)){
      if(newVal>pawianParams->UpperLimit(i)) newVal=pawianParams->UpperLimit(i);
      else if(newVal<pawianParams->LowerLimit(i)) newVal=pawianParams->LowerLimit(i);
    }
    pawianParams->SetValue(i,newVal);
  }
}


void AdamMinimizer::printFitResult(double evtWeightSumData){
  InfoMsg << "\n\n********************** Final fit parameters *************************" << endmsg;  
  _bestPawianParams->print(std::cout, true);
  InfoMsg << "\nFinal LH: " << _bestLH;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // calculate AIC, BIC criteria and output selected wave contrib
  ///////////////////////////////////////////////////////////////////////////////////////////////
  unsigned int noOfFreeFitParams=_bestPawianParams->VariableParameters();
  
  double BICcriterion=2.*_bestLH+noOfFreeFitParams*log(evtWeightSumData);
  double AICcriterion=2.*_bestLH+2.*noOfFreeFitParams;
  double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)
    / (evtWeightSumData-noOfFreeFitParams-1);
  InfoMsg << "\nnoOfFreeFitParams:\t" <<noOfFreeFitParams << endmsg;
  InfoMsg << "evtWeightSumData:\t" <<evtWeightSumData << endmsg;
  InfoMsg << "BIC:\t" << BICcriterion << endmsg;
  InfoMsg << "AIC:\t" << AICcriterion << endmsg;
  InfoMsg << "AICc:\t" << AICccriterion << endmsg;
}

