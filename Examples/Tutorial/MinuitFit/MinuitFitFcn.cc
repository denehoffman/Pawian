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

#include "Examples/Tutorial/MinuitFit/MinuitFitFcn.hh"
#include "Examples/Tutorial/MinuitFit/MinuitFit.hh"
#include "ErrLogger/ErrLogger.hh"
#include <cassert>

using namespace ROOT::Minuit2;

MinuitFitFcn::MinuitFitFcn(std::shared_ptr<MinuitFit> minuitFit) :
  _minFitPtr(minuitFit)
{
  if (0==_minFitPtr) {
    Alert << "MinuitFit pointer is 0 !!!!" << endmsg; 
    exit(1);
  }
}

MinuitFitFcn::~MinuitFitFcn()
{
}

double MinuitFitFcn::operator()(const std::vector<double>& par) const
{
  double result=_minFitPtr->calcChiSqr(par);
  DebugMsg << "current chi^2:\t"<< result << endmsg;  
  return result;
}

double MinuitFitFcn::Up() const 
{
return 1.;
}



