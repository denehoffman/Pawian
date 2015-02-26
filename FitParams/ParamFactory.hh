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

// ParamFactory class definition file. -*- C++ -*-
// Copyright 2015 Bertram Kopf

#pragma once

#include <iostream>
// #include <vector>
// #include <complex>
// #include <map>
#include <string>

// #include <cassert>
#include <memory>
#include "FitParams/AbsPawianParameters.hh"

//class AbsPawianParameters;

class ParamFactory{

public:

  /** Destructor */
  virtual ~ParamFactory();

  static ParamFactory* instance();
  std::shared_ptr<AbsPawianParameters> getParametersPointer(std::string type);

protected:
 ///Constructor 
  ParamFactory();
  static ParamFactory* _instance;

private:
};



