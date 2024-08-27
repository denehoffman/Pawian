//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

#include <math.h>
#include <stdio.h>
#include <boost/timer/timer.hpp>

#include "MinFunctions/AbsFcn.hh"
#include "Minuit2/FCNBase.h"
#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/MnUserParameters.h"

#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/NetworkServer.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

boost::timer::cpu_timer theTimer1;
boost::timer::cpu_timer theTimerAll;

template<typename T>
AbsFcn<T>::AbsFcn() :
  _fcnCounter(0)
  , _currentResFileName("currentResult"+GlobalEnv::instance()->outputFileNameSuffix()+".dat")
  ,_currentLHStreamFileName("currentLHVals"+GlobalEnv::instance()->outputFileNameSuffix()+".dat")
  ,_LHStream(new std::ofstream(_currentLHStreamFileName, ios::out))
{
}

template<typename T>
AbsFcn<T>::~AbsFcn()
{
}

template<typename T>
double AbsFcn<T>::Up() const
{
return .5;
}

template<typename T>
void AbsFcn<T>::printTimer() const{

  theTimer1.stop();
  boost::timer::cpu_times elapsed(theTimer1.elapsed());
  if(elapsed.wall > 0){
    InfoMsg << "Wall time: " << elapsed.wall / 1E9 << "s User: "
	    << elapsed.user/1E9 << "s System: " << elapsed.system/1E9 <<"s\n" << endmsg;
  }
  boost::timer::cpu_times elapsedAll(theTimerAll.elapsed());
  if(elapsedAll.wall > 0){
    InfoMsg << "Wtotal time: " << elapsedAll.wall / 1E9 << "s User: "
	    << elapsedAll.user/1E9 << "s System: " << elapsedAll.system/1E9 << "s\n" << endmsg;
  }
  theTimer1.start();
}

template<typename T>
void AbsFcn<T>::printFitParams(std::shared_ptr<AbsPawianParameters> par) const{
    par->print(std::cout, true);
}

template<typename T>
void  AbsFcn<T>::dumpFitParams(std::shared_ptr<AbsPawianParameters> par) const{
    std::ofstream theStream (_currentResFileName.c_str());
    par->print(theStream);
}

template<typename T>
void AbsFcn<T>::dumpLhVals(std::string input) const{
  if (!_LHStream->is_open()) _LHStream->open(_currentLHStreamFileName, std::ios_base::app);
  (*_LHStream) << _fcnCounter << "\t" << input << "\n";
  _LHStream->close();
}

template AbsFcn<FCNBase>::AbsFcn();
template AbsFcn<FCNBase>::~AbsFcn();
template AbsFcn<FCNGradientBase>::AbsFcn();
template AbsFcn<FCNGradientBase>::~AbsFcn();
