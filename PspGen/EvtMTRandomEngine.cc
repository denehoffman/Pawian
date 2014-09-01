//************************************************************************//
//                                                                        //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)                      //
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

// Copyright 2014 Julian Pychy

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "PspGen/EvtMTRandomEngine.hh"




double EvtMTRandomEngine::random(){
   return _variateGenerator();
}



EvtMTRandomEngine::EvtMTRandomEngine() :
   _seed(1),
   _engine(_seed),
   _variateGenerator(_engine, boost::random::uniform_real_distribution<>(0.,1.))
{
}



EvtMTRandomEngine::EvtMTRandomEngine(unsigned long int seed) :
   _seed(seed),
   _engine(_seed),
   _variateGenerator(_engine, boost::random::uniform_real_distribution<>(0.,1.))
{
}



void EvtMTRandomEngine::reset(){
   _variateGenerator.engine().seed(_seed);
   _variateGenerator.distribution().reset();
}
