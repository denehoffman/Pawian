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

#pragma once

// #include <iostream>
// #include <fstream>
// #include <string>
#include <vector>
// #include <map>
// #include <cassert>
// #include <memory>

// #include "TROOT.h"
// #include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class BBUnstableParFit {

public:

  BBUnstableParFit(double _mRes, double _wRes, double _m1, double _m2, std::string CMname="Dudek");
  virtual ~BBUnstableParFit();

  double calcNormD(const std::vector<double>& minPar);
  double calcNormDInvNeg(const std::vector<double>& minPar);
 
private:
  double mRes;
  double wRes;
  double m1;
  double m2;
  std::string _cmName;
};
