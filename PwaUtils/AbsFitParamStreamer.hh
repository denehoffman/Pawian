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

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/DataUtils.hh"

using namespace std;

typedef map<string, pair<double, double> > StringPairMap;

class AbsFitParamStreamer {
public:
  AbsFitParamStreamer(std::string& filePath);
  virtual ~AbsFitParamStreamer();

protected:
  StringPairMap _stringPairMap;
  virtual void fillParamMap() = 0;
  virtual void StringParameterMap(const string &theName, const double firstVal,
				  const double secondVal);
  virtual void fillParamMapAmps(std::vector< std::shared_ptr<const JPCLS> >& theJPCLSs, 
				std::string& suffix, std::map< std::shared_ptr<const JPCLS>,
				pair<double, double>, 
				pawian::Collection::SharedPtrLess >& toFill);
  virtual void fillParamMapMass(std::string& name, pair<double, double>& toFill);
  virtual void fillParamFlatte(std::string& name, double& toFill);
};
