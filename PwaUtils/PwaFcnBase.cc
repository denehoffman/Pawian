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

#include "Minuit2/MnUserParameters.h"

#include "PwaUtils/PwaFcnBase.hh"
#include "PwaUtils/AbsLh.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;


PwaFcnBase::PwaFcnBase(std::shared_ptr<AbsLh> absLh, std::shared_ptr<FitParamsBase> fitParamsBase, std::string suffix) :
  AbsFcn(fitParamsBase, suffix)
  ,_absLhPtr(absLh)
{
   if (0==_absLhPtr) { Alert << "AbsLh* _absLhPtr pointer is 0 !!!!" << endmsg; exit(1); }
   _absLhPtr->getDefaultParams(_defaultFitValParms, _defaultFitErrParms);
  
}

PwaFcnBase::~PwaFcnBase()
{
}

double PwaFcnBase::operator()(const std::vector<double>& par) const
{
  double result=0;

  fitParams theFitParmValTmp=_defaultFitValParms;

  _fitParamsBasePtr->getFitParamVal(par, theFitParmValTmp);
    
  result = _absLhPtr->calcLogLh(theFitParmValTmp);


  _fcnCounter++;

  if(_fcnCounter%20 == 0) printTimer();  
  printFitParams(par);
  dumpFitParams(par);

  return result;
}




