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

// AbsParamHandler class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <float.h>
#include "FitParams/AbsParamHandler.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "ErrLogger/ErrLogger.hh"

AbsParamHandler::AbsParamHandler() :
  _cacheAmps(false)
  ,_recalculate(true)
{
}

AbsParamHandler::~AbsParamHandler()
{
}


bool AbsParamHandler::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, 
					 std::shared_ptr<AbsPawianParameters> fitParOld){
  _recalculate=false;
  if(!_cacheAmps) _recalculate=true;
  else{ 
    std::vector<std::string>::iterator it;
    for(it=_paramNameList.begin(); it!=_paramNameList.end(); ++it){
      if (!CheckDoubleEquality(fitParNew->Value(*it), fitParOld->Value(*it))){
	_recalculate=true;
	return _recalculate;
      }
    } 
  }
  return _recalculate;
}

bool AbsParamHandler::CheckDoubleEquality(double a, double b){

   float diff = fabs(a - b);
   if (diff <= 10*DBL_EPSILON)
      return true;

   a = fabs(a);
   b = fabs(b);
   float largest = (b > a) ? b : a;

   if (diff <= largest * 2*DBL_EPSILON)
      return true;
   return false;
}
