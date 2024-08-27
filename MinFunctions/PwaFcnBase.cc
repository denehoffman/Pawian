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
#include <iomanip>
#include <boost/timer/timer.hpp>

#include "Minuit2/FCNBase.h"
#include "MinFunctions/PwaFcnBase.hh"
#include "FitParams/ParamDepHandler.hh"

#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ConfigParser/ParserBase.hh"
#include "ErrLogger/ErrLogger.hh"

template<typename T>
PwaFcnBase<T>::PwaFcnBase() :
  AbsFcn<T>()
{
   this->_currentPawianParms = GlobalEnv::instance()->startPawianParams();
}

template<typename T>
PwaFcnBase<T>::~PwaFcnBase()
{
}

template<typename T>
double PwaFcnBase<T>::operator()(const std::vector<double>& par) const
{
  double result=0;
  this->_currentPawianParms->SetAllValues(par);
  ParamDepHandler::instance()->ApplyDependencies(this->_currentPawianParms);

  result = GlobalEnv::instance()->Channel()->Lh()->calcLogLh(this->_currentPawianParms);

  if(this->_fcnCounter%GlobalEnv::instance()->parser()->stepSizeLhPrint() == 0){
    InfoMsg << "current LH = " << std::setprecision(16) << result << endmsg;
  }
  if(this->_fcnCounter%GlobalEnv::instance()->parser()->stepSizeTimer() == 0) this->printTimer();
  if(this->_fcnCounter%GlobalEnv::instance()->parser()->stepSizeParamsPrint() == 0) this->printFitParams(this->_currentPawianParms);
  if(this->_fcnCounter%GlobalEnv::instance()->parser()->stepSizeParamsDump() == 0) this->dumpFitParams(this->_currentPawianParms);
  if(this->_fcnCounter%GlobalEnv::instance()->parser()->stepSizeLhDump() == 0){
    std::string resultString = std::to_string(result);
    this->dumpLhVals(resultString);
  }
  this->_fcnCounter++;
  return result;
}

template PwaFcnBase<FCNBase>::PwaFcnBase();
template PwaFcnBase<FCNBase>::~PwaFcnBase();




