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
   Info << evologo << endmsg;
}


std::vector<double> EvoMinimizer::Minimize(){
  
   double startlh = (*_theAbsFcn)(_bestParamsGlobal.Params());
   double minlh = startlh;
   double itlh = startlh;
   int numnoimprovement = 0;

   Info << "Start LH = " << startlh << endmsg;

   for(int i = 0; i<_iterations; i++){

      int numbetterlh = 0;
      double maxitlhspread=0;

      for(int j = 0; j<_population; j++){
	 _tmpParams = _bestParamsGlobal;
	 itlh = minlh;
	 ShuffleParams();
	 double currentlh = (*_theAbsFcn)(_tmpParams.Params());
	 if(fabs(currentlh - minlh) > maxitlhspread){
	    maxitlhspread = fabs(currentlh - minlh);
	 }
	 if(currentlh < minlh){
	    minlh = currentlh;
	    _bestParamsGlobal = _tmpParams;
	 }
	 if(currentlh < itlh){
	    numbetterlh++;
	 }	 
      }

      if(maxitlhspread < 0.1)
	 break;

      if(numbetterlh == 0)
	 numnoimprovement++;

      if(numnoimprovement>5){
	 numnoimprovement=0;
	 AdjustSigma(0.5);
      }
      
      Info << "===============================================" << endmsg;
      Info << "Iteration " << i+1 << " / " << _iterations << " done. Best LH " << minlh << endmsg;
      Info << "Likelihood improvements " << numbetterlh << " / " << _population << endmsg;

      if(((double)numbetterlh / (double)_population) <= 0.01){
	 AdjustSigma(0.9);
	 Info << "Decreasing errors." << endmsg;
      }
      else if(((double)numbetterlh / (double)_population) > 0.05){
	 AdjustSigma(1.1);
	 Info << "Increasing errors." << endmsg;
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



void EvoMinimizer::AdjustSigma(double factor){

   for(unsigned int i=0; i<_bestParamsGlobal.Params().size(); i++){
      if(_bestParamsGlobal.Parameter(i).IsFixed())
	 continue;

      _bestParamsGlobal.SetError(i, _bestParamsGlobal.Error(i) * factor);
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
