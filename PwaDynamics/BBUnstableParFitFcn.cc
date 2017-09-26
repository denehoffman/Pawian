//************************************************************************//
//									  //
//  Copyright 2017 Markus Kuhlmann (mkuhlmann@ep1.rub.de)		  //
//  	      	   - Ruhr-Universität Bochum 				  //
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

#include "PwaDynamics/BBUnstableParFitFcn.hh"
#include "PwaDynamics/BBUnstableParFit.hh"
#include "ErrLogger/ErrLogger.hh"
#include <cassert>

using namespace ROOT::Minuit2;

BBUnstableParFitFcn::BBUnstableParFitFcn(std::shared_ptr<BBUnstableParFit> minuitFit) :
  _minFitPtr(minuitFit)
{
  if (0==_minFitPtr) {
    Alert << "BBUnstableParFit pointer is 0 !!!!" << endmsg; 
    exit(1);
  }
}

BBUnstableParFitFcn::~BBUnstableParFitFcn()
{
}

double BBUnstableParFitFcn::operator()(const std::vector<double>& par) const
{
  double result=_minFitPtr->calcNormD(par);
  DebugMsg << "current value:\t"<< result << endmsg;  
  return result;
}

double BBUnstableParFitFcn::Up() const 
{
return 0.5;
}



