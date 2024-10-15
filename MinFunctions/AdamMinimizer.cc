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

AdamMinimizer::AdamMinimizer(std::shared_ptr<AbsFcn<FCNGradientBase>> theAbsFcnPtr, std::shared_ptr<AbsPawianParameters> upar)
    : AbsPawianMinimizer<FCNGradientBase>(theAbsFcnPtr, upar),
      _max_iterations(5000),
      _s(upar->Params().size(), 0.0),
      _v(upar->Params().size(), 0.0),
      _learning_rate(0.001),
      _initial_lr(0.001),
      _decay_rate(0.01),
      _iterations(0),
      _currentPawianParams(std::shared_ptr<AbsPawianParameters>(upar->Clone())),
      _bestLH(1.e20),
      _bestPawianParams(std::shared_ptr<AbsPawianParameters>(upar->Clone())),
      _noItersWoImprovement(0),
      _patience(10) { 
}

AdamMinimizer::~AdamMinimizer() {}

void AdamMinimizer::minimize() {
 
    const double gradient_tolerance = 1e-8;

    if (_iterations == 0) {
        double currentLH = _absFcn->operator()(_bestPawianParams->Params());
        if (currentLH < _bestLH) _bestLH = currentLH;
        std::cout << "Initial best LH: " << _bestLH << "\tnoIters: " << _iterations << "\tnoItersWoImprovement: " << _noItersWoImprovement << std::endl;
    }

    while (_iterations <= _max_iterations && _noItersWoImprovement < 1000) {
        if (_noItersWoImprovement > 0 && _noItersWoImprovement >= _patience) {
            _currentPawianParams = std::shared_ptr<AbsPawianParameters>(_bestPawianParams->Clone());
	    _patience = std::min(_patience * 1.05, 30.0);
	    _noItersWoImprovement = 0;
        }
        std::vector<double> derivatives = _absFcn->Gradient(_currentPawianParams->Params());

        updateParameters(_currentPawianParams, derivatives, _s, _v, _iterations);

        double currentLH = _absFcn->operator()(_currentPawianParams->Params());

        if (currentLH < _bestLH) {
            _bestLH = currentLH;
            _bestPawianParams = std::shared_ptr<AbsPawianParameters>(_currentPawianParams->Clone());
            _noItersWoImprovement = 0;
        } else {
            ++_noItersWoImprovement;
        }

        double gradientNorm = 0.0;
        for (const auto& grad : derivatives) {
            gradientNorm += grad * grad;
        }
        gradientNorm = sqrt(gradientNorm);
	
	if (gradientNorm < gradient_tolerance && std::abs(currentLH - _bestLH) / std::abs(_bestLH) < 1e-4) {
	  std::cout << "Stopping early: Gradient norm below tolerance and no significant likelihood improvement." << std::endl;
	  break;
	}


        ++_iterations;
	_learning_rate *= 0.999;
	std::cout << "best LH: " << _bestLH << "\tnoIters: " << _iterations << "\tnoItersWoImprovement: " << _noItersWoImprovement << std::endl;
    }
}

void AdamMinimizer::updateParameters(std::shared_ptr<AbsPawianParameters> pawianParams, std::vector<double>& gradients, std::vector<double>& s, std::vector<double>& v, int t) {
    double beta1 = 0.9;
    double beta2 = 0.999;
    double epsilon = 1e-8;

    for (unsigned int i = 0; i < pawianParams->Params().size(); ++i) {
        if (pawianParams->IsFixed(i)) continue;
            
	double clip_value = 0.2 * std::abs(pawianParams->Value(i));
        gradients[i] = std::min(clip_value, std::max(-clip_value, gradients[i]));

        s.at(i) = beta1 * s.at(i) + (1.0 - beta1) * gradients.at(i);
        v.at(i) = beta2 * v.at(i) + (1.0 - beta2) * gradients.at(i) * gradients.at(i);

        
        double s_hat = s.at(i) / (1.0 - pow(beta1, (t + 1)));
        double v_hat = v.at(i) / (1.0 - pow(beta2, (t + 1)));


	//double learning_rate_t = _initial_lr / (1.0 + _decay_rate * t);
	//double adaptive_lr = learning_rate_t / (std::abs(pawianParams->Value(i)) + epsilon);
	//double gradient_update = adaptive_lr * s_hat / (std::sqrt(v_hat) + epsilon);
	double newVal = pawianParams->Value(i) - _learning_rate * s_hat / (std::sqrt(v_hat) + epsilon);
       
	
        if (pawianParams->HasLimits(i)) {
            if (newVal > pawianParams->UpperLimit(i)) newVal = pawianParams->UpperLimit(i);
            else if (newVal < pawianParams->LowerLimit(i)) newVal = pawianParams->LowerLimit(i);
        }

	pawianParams->SetValue(i, newVal);
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

