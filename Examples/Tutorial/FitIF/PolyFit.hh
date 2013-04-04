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

#ifndef _PolyFit_H
#define _PolyFit_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "Examples/Tutorial/FitIF/MinimizerInterface/MIData.hh"

class TFile;
class TGraph;
class TCanvas;
class TRandom;

class PolyFit : public MIData {

public:

  // create/copy/destroy:

  ///Constructor 
  PolyFit(double p0, double p1, double p2, double p3, double sigma);


  /** Destructor */
  virtual ~PolyFit();


  double controlParameter(const std::vector<double>& minPar);
  void drawGraph(double a, double b, double c, double d);
  // Getters:
 
protected:


private:
  TFile* _theTFile;
  std::map <unsigned int, TGraph* > _myGraph;

  std::vector< double > _xValue;
  std::vector< double > _yValue;

  double _sigma;

};

#endif
