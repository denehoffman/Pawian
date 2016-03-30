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

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

class AbsPawianParameters;

class AbsParamHandler {

public:
  AbsParamHandler();
  virtual ~AbsParamHandler();

  // virtual void getDefaultParams(fitParCol& fitVal, fitParCol& fitErr)=0;
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar)=0;
  // virtual bool checkRecalculation(fitParCol& theParamVal)=0;
  virtual bool checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld);
  virtual void cacheAmplitudes()=0;
  // virtual void updateFitParams(fitParCol& theParamVal)=0;
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){;}
  static bool CheckDoubleEquality(double a, double b);
  std::vector<std::string> paramNames() {return _paramNameList;}
  virtual void fillParamNameList(){;}
protected:
  bool _cacheAmps;
  bool _recalculate;
  std::vector<std::string> _paramNameList;

private:

};
