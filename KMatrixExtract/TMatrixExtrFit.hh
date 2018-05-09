//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
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
//TMatrixExtrFit class definition file. -*- C++ -*-
// Copyright 2018 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <map>
#include <memory>

#include "KMatrixExtract/TMatrixExtrBase.hh"

class TMatrixExtrFit : public TMatrixExtrBase {

public:

  // create/copy/destroy:

  ///Constructor 
  TMatrixExtrFit(std::string pathToConfigParser, std::string pathToFitParams, std::complex<double> energyBorderMin, std::complex<double> energyBorderMax, std::complex<double> energyStartParams);


  /** Destructor */
  virtual ~TMatrixExtrFit();

  // Getters:
  double calcTMatrix(double eReal, double eImag);
  std::complex<double> energyMin() {return _energyMin;}
  std::complex<double> energyMax() {return _energyMax;}
  std::complex<double> energyStart() {return _energyStart;}

protected:


private:
  std::complex<double> _energyMin;
  std::complex<double> _energyMax;
  std::complex<double> _energyStart;
};


