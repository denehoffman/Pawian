//************************************************************************//
//                                                                        //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 Julian Pychy (julian@ep1.rub.de)                       //
//                 - Ruhr-Universit??t Bochum                              //
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

#define mysign(x) (( x > 0 ) - ( x < 0 ))

#include <iostream>
#include <fstream>
#include <boost/random.hpp>

#include "MinFunctions/EvoMinimizer.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ConfigParser/ParserBase.hh"

const double EvoMinimizer::DECREASESIGMAFACTOR = 0.9;
const double EvoMinimizer::INCREASESIGMAFACTOR = 1.1;
const double EvoMinimizer::DECREASELOWTHRESH = 0.2;
const double EvoMinimizer::INCREASEHIGHTHRESH = 0.2;
const double EvoMinimizer::LHSPREADEXIT = 0.01;


// Constructor takes AbsFcn to minimze, start parameters upar, population and iteration
// sizes and the output file name suffix
EvoMinimizer::EvoMinimizer(std::shared_ptr<AbsFcn> theAbsFcnPtr, 
			   std::shared_ptr<AbsPawianParameters> upar, 
			   int population, int iterations) :
  AbsPawianMinimizer(theAbsFcnPtr, upar)
  ,_population(population)
  , _iterations(iterations)
  , _evoRatioOfModParams(GlobalEnv::instance()->parser()->evoRatioOfModParams())
  , _currentResultFileName("currentEvoResult"+GlobalEnv::instance()->outputFileNameSuffix()+".dat")
{
  if (_evoRatioOfModParams <= 0. || _evoRatioOfModParams>1.){
    Alert << "_evoRatioOfModParams = " << _evoRatioOfModParams << " not possible\n"
  	  << "value must be set between 0. and 1. !!!!" << endmsg;
    exit(1);
  }
  InfoMsg << "population: " << _population
       <<"\niterations: " << _iterations
       << "\nratio of parameters to be modified: " << _evoRatioOfModParams << endmsg;
}


// Minimization takes place here
void EvoMinimizer::minimize(){

   double startlh = (*_absFcn)(_bestPawianParams->Params());
   double minlh = startlh;
   int numnoimprovement = 0;

   InfoMsg << "Start LH = " << startlh << endmsg;

   for(int i = 0; i<_iterations; i++){

      int numbetterlh = 0;
      double maxitlhspread=0;
      double itlh = minlh;
      _iterationParamBackup = std::shared_ptr<AbsPawianParameters>(_bestPawianParams->Clone());
      _bestParamsIteration = std::shared_ptr<AbsPawianParameters>(_bestPawianParams->Clone());

      for(int j = 0; j<_population; j++){

         // Get iteration start parameters and shuffle them
	_tmpParams = std::shared_ptr<AbsPawianParameters>(_bestPawianParams->Clone());
         ShuffleParams();

         // Calc likelihood
         double currentlh = (*_absFcn)(_tmpParams->Params());

         // Get information for break condition
         if(fabs(currentlh - minlh) > maxitlhspread){
               maxitlhspread = fabs(currentlh - minlh);
         }

         // Test for new best lh in iteration;
         if(currentlh < itlh){
            itlh = currentlh;
            _bestParamsIteration = std::shared_ptr<AbsPawianParameters>(_tmpParams->Clone());
         }

         // Count number of lh improvements in the iteration
         if(currentlh < minlh){
            numbetterlh++;
        }
      }

      // Break condition
      if(maxitlhspread < LHSPREADEXIT)
         break;

      // If a new minimum has been found, store information
      // and print parameters
      if(numbetterlh > 0){
	_bestPawianParams = std::shared_ptr<AbsPawianParameters>(_bestParamsIteration->Clone());
         minlh = itlh;
         numnoimprovement=0;

         std::ofstream theStream(_currentResultFileName.c_str());
	_bestPawianParams->print(theStream);
       }
      else{
         numnoimprovement++;
      }

      // Emergency sigma decrease
      if(numnoimprovement > 5){
         AdjustSigma(DECREASESIGMAFACTOR*DECREASESIGMAFACTOR, 0);
      }

      // Adjust errors depending on ratio of improved likelihoods to population size
      std::string logmessage;
      if(((double)numbetterlh / (double)_population) < DECREASELOWTHRESH){
         AdjustSigma(DECREASESIGMAFACTOR, numbetterlh);
         logmessage = "Decreasing errors.";
      }
      else if(((double)numbetterlh / (double)_population) >= INCREASEHIGHTHRESH){
         AdjustSigma(INCREASESIGMAFACTOR, numbetterlh);
         logmessage = "Increasing errors.";
      }

      // Print iteration summary
      InfoMsg << "===============================================" << endmsg;
      InfoMsg << "Iteration " << i+1 << " / " << _iterations << " done. Best LH " << minlh << endmsg;
      InfoMsg << "Likelihood improvements " << numbetterlh << " / " << _population << endmsg;
      InfoMsg << "Likelihood spread " << maxitlhspread << endmsg;
      InfoMsg << logmessage << endmsg;
      InfoMsg << "===============================================" << endmsg;

   } // Iterations

   //   return _bestParamsGlobal->Params();
   _minimumReached=true;
}




// Shuffles parameters using gauss distributions
void EvoMinimizer::ShuffleParams(){

   typedef boost::normal_distribution<double> NormalDistribution;
   typedef boost::mt19937 RandomGenerator;
   typedef boost::variate_generator<RandomGenerator&,	NormalDistribution> GaussianGenerator;
   static RandomGenerator rng(static_cast<unsigned> (time(0)));

   bool acceptNewParams=false;
   while(!acceptNewParams){
     for(unsigned int i=0; i<_tmpParams->Params().size(); i++){
       
       // Don't touch fixed parameters
       if(_tmpParams->IsFixed(i)){
         continue;
       }

       boost::random::uniform_real_distribution<> coinReal(0.,1.);
       double c01 = coinReal(rng);
       if (c01 > _evoRatioOfModParams) continue; //only 10% of the parameters changed 
       acceptNewParams=true;
       
       // Decide whether parameter is increased or decreased
       boost::random::uniform_int_distribution<> coin(0,1);
       bool c = coin(rng);
       
       // Initialize gaussian width as parameter error
       double sigma = _tmpParams->Error(i);
       
       // If gaussian collides with parameter limit, reduce width
       if(_tmpParams->HasLimits(i)){
         if(c && (2*sigma > (_tmpParams->UpperLimit(i) - _tmpParams->Value(i)))){
	   sigma = (_tmpParams->UpperLimit(i) - _tmpParams->Value(i)) / 2.0;
         }
         else if(!c && (2*sigma > (_tmpParams->Value(i) - _tmpParams->LowerLimit(i)))){
	   sigma = (_tmpParams->Value(i) - _tmpParams->LowerLimit(i)) / 2.0;
         }
       }
       
       // Get random number and set new parameter
       NormalDistribution gaussian_dist(0, sigma);
       GaussianGenerator generator(rng, gaussian_dist);
       double val = fabs(generator());
       
       if(c) _tmpParams->SetValue(i, _tmpParams->Value(i) + val);
       else  _tmpParams->SetValue(i, _tmpParams->Value(i) - val);
       
       // Check for limits
       if(_tmpParams->HasLimits(i)){
	 if(_tmpParams->Value(i) < _tmpParams->LowerLimit(i))
	   _tmpParams->SetValue(i, _tmpParams->LowerLimit(i));
	 if(_tmpParams->Value(i) > _tmpParams->UpperLimit(i))
	   _tmpParams->SetValue(i, _tmpParams->UpperLimit(i));
       }
     }
   }
}

void EvoMinimizer::printFitResult(double evtWeightSumData){
  if(!_minimumReached){
    Alert << "minimum has not been reached!!!" << endmsg;
    exit(1);
  }

  InfoMsg <<"Final fit params:" << endmsg;
  std::vector<double> pVals=_bestPawianParams->Params();  
  std::vector<double> pErrs=_bestPawianParams->Errors();
  for(unsigned int i=0; i < pVals.size(); ++i){
    InfoMsg << _bestPawianParams->GetName(i) << " : " << pVals.at(i) << "\t" << pErrs.at(i) << endmsg;
  }

  double finalLh = (*_absFcn)(_bestPawianParams->Params());
  InfoMsg <<"final result theLh = "<< finalLh << endmsg;
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // calculate AIC, BIC criteria and output selected wave contrib
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  unsigned int noOfFreeFitParams=_bestPawianParams->VariableParameters();

  double BICcriterion=2.*finalLh+noOfFreeFitParams*log(evtWeightSumData);
  double AICcriterion=2.*finalLh+2.*noOfFreeFitParams;
  double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/(evtWeightSumData-noOfFreeFitParams-1);
  InfoMsg << "noOfFreeFitParams:\t" << noOfFreeFitParams << endmsg;
  InfoMsg << "evtWeightSumData:\t" << evtWeightSumData << endmsg;
  InfoMsg << "BIC:\t" << BICcriterion << endmsg;
  InfoMsg << "AIC:\t" << AICcriterion << endmsg;
  InfoMsg << "AICc:\t" << AICccriterion << endmsg;
}



// Increase or decrease parameter errors by a factor
void EvoMinimizer::AdjustSigma(double factor, int numimprovements){

   for(unsigned int i=0; i<_bestPawianParams->Params().size(); i++){

      // Don't touch fixed parameters
      if(_bestPawianParams->IsFixed(i)){
         continue;
      }

      // When a lh improvement was achieved, don't decrease errors of parameters
      // that changed rapidly, also add a bonus to increasements.
      double afactor=factor;
      double pardiffsigmas = fabs(_bestPawianParams->Value(i) - _iterationParamBackup->Value(i)) / _bestPawianParams->Error(i);

      if(numimprovements > 0 && pardiffsigmas > 1.5){
         if(factor < 1)        afactor = 1.0;
         else if(factor > 1 )  afactor *= afactor;
      }

      // Set new error
      _bestPawianParams->SetError(i, _bestPawianParams->Error(i) * afactor);
   }

}

