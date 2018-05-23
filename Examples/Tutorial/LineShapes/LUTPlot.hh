//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//		   Markus Kuhlmann (mkuhlmann@ep1.rub.de)		  //
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
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cassert>

#include <memory>

#include "TROOT.h"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class TFile;
class TH1F;
class TH2F;

class LUTPlot {

public:

  // create/copy/destroy:

  ///Constructor 
  LUTPlot(std::string LUTFilePath, double realGranularity, double imagGranularity);


  /** Destructor */
  virtual ~LUTPlot();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  
  TH2F* _RealPart2D;
  TH2F* _ImagPart2D;

  TH1F* _RealPart;
  TH1F* _ImagPart;
  
  const std::string _LUTFilePath;
  const double _realGranularity;
  const double _imagGranularity;
};


