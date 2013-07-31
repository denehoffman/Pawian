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

#define mysign(x) (( x > 0 ) - ( x < 0 ))

#include <iostream>
#include <fstream>
#include <boost/random.hpp> 
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/EvoMinimizer.hh"
#include "PwaUtils/FitParamsBase.hh"

const double EvoMinimizer::DECREASESIGMAFACTOR = -0.1;
const double EvoMinimizer::INCREASESIGMAFACTOR =  0.1;
const double EvoMinimizer::DECREASELOWTHRESH = 0.001;
const double EvoMinimizer::INCREASEHIGHTHRESH = 0.05;
const double EvoMinimizer::LHSPREADEXIT = 0.1;



EvoMinimizer::EvoMinimizer(AbsFcn& theAbsFcn, MnUserParameters upar, int population, int iterations, std::string suffix) :
  _population(population)
  , _iterations(iterations)
  , _theAbsFcn(&theAbsFcn)
  , _fitParamBase(theAbsFcn.fitParamBase())
  , _currentBestParams(theAbsFcn.defaultFitValParms())
  , _defaultFitErrParms(theAbsFcn.defaultFitErrParms())
  , _currentResultFileName("currentEvoResult"+suffix+".dat")
{
   _bestParamsGlobal = upar;
   _parShuffleMod.assign(upar.Parameters().size(), 0);
   Info << evologo << endmsg;
}



std::vector<double> EvoMinimizer::Minimize(){
  
   double startlh = (*_theAbsFcn)(_bestParamsGlobal.Params());
   double minlh = startlh;
   int numnoimprovement = 0;

   Info << "Start LH = " << startlh << endmsg;

   for(int i = 0; i<_iterations; i++){

      int numbetterlh = 0;
      double maxitlhspread=0;
      double itlh = minlh;

      _iterationParamBackup = _bestParamsGlobal;
      _bestParamsIteration = _bestParamsGlobal;

      for(int j = 0; j<_population; j++){
	 _tmpParams = _bestParamsGlobal;
	 ShuffleParams();
	 double currentlh = (*_theAbsFcn)(_tmpParams.Params());
	 if(fabs(currentlh - minlh) > maxitlhspread){
	    maxitlhspread = fabs(currentlh - minlh);
	 }
	 if(currentlh < itlh){
	    itlh = currentlh;
	    _bestParamsIteration = _tmpParams;
	    _fitParamBase->getFitParamVal(_bestParamsGlobal.Params(), _currentBestParams); 
	    std::ofstream theStream(_currentResultFileName.c_str());
	    _fitParamBase->dumpParams(theStream, _currentBestParams, _defaultFitErrParms); 
	 }
	 if(currentlh < minlh){
	    numbetterlh++;
	 }	 
      }


      if(maxitlhspread < LHSPREADEXIT)
	 break;

      if(numbetterlh > 0){
	 _bestParamsGlobal = _bestParamsIteration;
	 minlh = itlh;
	 numnoimprovement=0;
      }
      else
	 numnoimprovement++;


      if(numnoimprovement>5){
	 AdjustSigma(DECREASESIGMAFACTOR*2, 0);
      }
      
      Info << "===============================================" << endmsg;
      Info << "Iteration " << i+1 << " / " << _iterations << " done. Best LH " << minlh << endmsg;
      Info << "Likelihood improvements " << numbetterlh << " / " << _population << endmsg;

      if(((double)numbetterlh / (double)_population) <= DECREASELOWTHRESH){
	 AdjustSigma(DECREASESIGMAFACTOR, numbetterlh);
	 Info << "Decreasing errors." << endmsg;
      }
      else if(((double)numbetterlh / (double)_population) > INCREASEHIGHTHRESH){
	 AdjustSigma(INCREASESIGMAFACTOR, numbetterlh);
	 Info << "Increasing errors." << endmsg;
      }
      else{
	 AdjustSigma(0, numbetterlh);
      }
      Info << "===============================================" << endmsg;
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



void EvoMinimizer::AdjustSigma(double factor, int numimprovements){

   for(unsigned int i=0; i<_bestParamsGlobal.Params().size(); i++){
      if(_bestParamsGlobal.Parameter(i).IsFixed())
	 continue;

      if(numimprovements>0){
	 double pardiff = fabs(_bestParamsGlobal.Value(i) - _iterationParamBackup.Value(i))
	    / _bestParamsGlobal.Error(i);

	 if(pardiff > 1.0 && _parShuffleMod[i] < 1.0)
	    _parShuffleMod[i] += 0.05; 
	 else if(pardiff < 0.4 && _parShuffleMod[i] > -1.0)
	    _parShuffleMod[i] -= 0.05; 
      }

      double mod = factor * pow( fabs(_parShuffleMod[i]) + 1, mysign(_parShuffleMod[i]));
      _bestParamsGlobal.SetError(i, _bestParamsGlobal.Error(i) * (1 + mod));
   }
}



const char* EvoMinimizer::evologo = R"(
                                   ....$$$$ZZZ$$Z7ZZ:..... .                                                            
                                   .~$ZZ.............O$ZZ+.                                                             
         .. ... ...=...............ZZZ.. .... ..      ..ZOO: ...                                                        
       ......$OZZ$7$ZZ$ZZZZZZZZZZ$ZZ.....ZO$....      .. .ZO. .                                                         
       ..:8OO...................:ZZ.......Z$. ..     ..  ..OO...                                                        
    ...OOO. ....                 ..     ....           ....$Z...                                                        
   .8OO..... .                          . ..            . +$Z...                                                        
. .OO...                                            .......$Z...                                                        
.OO. ...                                             ....Z:Z$... .                                                      
OO......  .     ........   ..                      ...$$:ZZ.ZZOZ.. .                                                    
OOO88O8OOI........... .. .......                .....ZZZZZ...$ZZO...                                                    
8O.....,$ZZOZZOOOO8OOOZ$........               ..$OZ$..........7$7..                                                    
88.. . .            ....~ZZOOZOOZZZ$....       .Z....ZZZ$...   .$$$:....                                                
.8O.....            ................ZZ..       .... ....ZZ..   ..$ZZ....                                                
..8OO...            .   ....... ........       ....... .OO$.    ..$:Z...                                                
 ...888~....                           .....       ....ZOZ..    ...ZZ$..                                                
......88OO. ..                         .....        .,OZZ..     ....ZZ$.                                                
       .888OO............. .         ..OOZ.. . ....ZZOO.           ..$$Z....                                            
        ...7888OO$........ .         ...ZOOZZZZDZZZO....           ...ZZ$...                                            
        .......ZO8O8OOZZ..          .  ...ZOOO8Z$..... .           ... ZZZ..                                            
                ....,ZOOOOZ..  .              ........=ZOOZZO$......   .ZZ$....                                        
                .........=O8....                 ..O88OZZ...$ZZZZ?..   ...ZZ....                                        
                       . ..88O... ..           ...OO?ZOZZ...  ..ZZZZ.......ZZ...                                        
                       ... .O8O?....            .O8,.ZOZZ$.   ....ZZZZ......ZO..                                        
                         .....ODOI..            .O8.....ZO..   .... .ZZZZ. ..OO.                                        
                            ..OOO8OI.. ..      ..O8.    ZZZ... .    ...Z$$$?.Z$....                                    
                            .O8..ZOOOOZ.....    ZZ..    .ZO....    .. ...ZZZ$$$$.. .                                    
                        ....O8,...ZOOO88O...... .Z..    ..ZZ....           .ZZ$$Z...
                        ...OO8.....ZOO.O8OZ. .....:.   .. .ZO...           ....$$Z$.
                       .. .Z8...... .OO$.ZOOO.Z.,O8.    . ..$O..            .. ..$$$.
                       ....OO........,OOOOOOOO$Z.OO.....   ..$O~.......         .ZZZ:
                      ....IZO.......8888O8OO..Z$ZZOO....    ..,Z8.......       ...ZO$O
               .......=O8888O.. .O88O....OO+OOOOZO88O........~,ZZZ$Z....... . ....ZZ..
                 .Z888OOO7....O88O............ZOZOOZOOOOZOOOO$.............. ....7ZOO...
              ..O888OOOO8OOOO8OZZ=IZZOZO8OOO8OOOZZOZ..ZZZ~$ZZ$?.7?I?::?$OOOZZZZZZ$Z.OOZ.
            . .O88888888OOOOOOOO8OOOOOOZZOOOZ~..,ZZ=.~=7......ZZZZ$$ZZZZ$Z$ZZZZZZZOZ... 
            .. ..........  .......    .Z$?ZZZOZ$IIZ$=7Z:=7ZZ~.....................$ZZO..      
                   . ...        ....   .......  .... ..$+I7Z$$ZZZZZZZ$ZZZZZZZZZZZZZ$Z$$..        
                                                      . ........ ...............+?=...        

//**********************************************************************************//
//                                                                                  //
//                           ApeFrog evolutionary minimizer                         //
//                                                                                  //
//**********************************************************************************//


 )";
