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

#include "MinFunctions/PwaFcnBase.hh"
#include "FitParams/ParamDepHandler.hh"

#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ConfigParser/ParserBase.hh"
#include "ErrLogger/ErrLogger.hh"

PwaFcnBase::PwaFcnBase() :
  AbsFcn()
{
   _currentPawianParms = GlobalEnv::instance()->startPawianParams();
}

PwaFcnBase::~PwaFcnBase()
{
}

double PwaFcnBase::operator()(const std::vector<double>& par) const
{
  double result=0;
  _currentPawianParms->SetAllValues(par);
  ParamDepHandler::instance()->ApplyDependencies(_currentPawianParms);

  result = GlobalEnv::instance()->Channel()->Lh()->calcLogLh(_currentPawianParms);

  if(_fcnCounter%GlobalEnv::instance()->parser()->stepSizeLhPrint() == 0){
    InfoMsg << "current LH = " << std::setprecision(16) << result << endmsg;
  }
  if(_fcnCounter%GlobalEnv::instance()->parser()->stepSizeTimer() == 0) printTimer();
  if(_fcnCounter%GlobalEnv::instance()->parser()->stepSizeParamsPrint() == 0) printFitParams(_currentPawianParms);
  if(_fcnCounter%GlobalEnv::instance()->parser()->stepSizeParamsDump() == 0) dumpFitParams(_currentPawianParms);
  if(_fcnCounter%GlobalEnv::instance()->parser()->stepSizeLhDump() == 0){
    std::string resultString = std::to_string(result);
    dumpLhVals(resultString);
  }
  _fcnCounter++;
  return result;
}




