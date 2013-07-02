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

#include <boost/random.hpp> 
#include "ErrLogger/ErrLogger.hh"
#include "EvoMinimizer.hh"

EvoMinimizer::EvoMinimizer(AbsFcn& theAbsFcn, MnUserParameters upar, int population, int iterations) :
   _population(population)
 , _iterations(iterations)
 , _theAbsFcn(&theAbsFcn)
{
   _bestParamsGlobal = upar;
}


std::vector<double> EvoMinimizer::Minimize(){
  
   double startlh = (*_theAbsFcn)(_bestParamsGlobal.Params());
   double minlh = startlh;
   double itlh = startlh;

   for(int i = 0; i<_iterations; i++){
      Info << "Iteration " << i+1 << " / " << _iterations << " Best LH "<< minlh << endmsg;
      _bestParamsIteration = _bestParamsGlobal;
      int numbetterlh = 0;

      for(int j = 0; j<_population; j++){
	 _tmpParams = _bestParamsIteration;
	 itlh = minlh;
	 ShuffleParams();
	 double currentlh = (*_theAbsFcn)(_tmpParams.Params());
	 if(currentlh < minlh){
	    minlh = currentlh;
	    _bestParamsGlobal = _tmpParams;
	 }
	 if(currentlh < itlh){
	    numbetterlh++;
	 }	 
      }
      
      if((numbetterlh / _population) < 0.15)
	 AdjustSigma(0.95);
      else if((numbetterlh / _population) > 0.2)
	 AdjustSigma(1.05);
   }
   
   return _bestParamsGlobal.Params();
}



void EvoMinimizer::ShuffleParams(){

   typedef boost::normal_distribution<double> NormalDistribution;
   typedef boost::mt19937 RandomGenerator;
   typedef boost::variate_generator<RandomGenerator&,		\
      NormalDistribution> GaussianGenerator;
   
   static RandomGenerator rng(static_cast<unsigned> (time(0)));
   

   for(unsigned int i=0; i<_tmpParams.Params().size(); i++){
      
      if(_tmpParams.Parameter(i).IsFixed()) continue;

      NormalDistribution gaussian_dist(_tmpParams.Parameter(i).Value(), _tmpParams.Parameter(i).Error());
      GaussianGenerator generator(rng, gaussian_dist);
      _tmpParams.SetValue(i,generator());

      if(_tmpParams.Parameter(i).HasLimits()){
	 if(_tmpParams.Value(i) < _tmpParams.Parameter(i).LowerLimit())
	    _tmpParams.SetValue(i, _tmpParams.Parameter(i).LowerLimit()); 
	 if(_tmpParams.Value(i) > _tmpParams.Parameter(i).UpperLimit())
	    _tmpParams.SetValue(i, _tmpParams.Parameter(i).UpperLimit()); 
       }
   }
}



void EvoMinimizer::AdjustSigma(double factor){

   for(unsigned int i=0; i<_bestParamsGlobal.Params().size(); i++){
      if(_bestParamsGlobal.Parameter(i).IsFixed())
	 continue;

      _bestParamsGlobal.SetError(i, _bestParamsGlobal.Error(i) * factor);
   }
}
