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

// PPole class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file PPole.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include <iostream>
#include <vector>
#include <complex>
#include <boost/shared_ptr.hpp>
using namespace std;

class AbsPhaseSpace;
//_____________________________________________________________________________
//_____________________________________________________________________________

class PPole : public KPole {

public:

  // create/copy/destroy:

  /// Default Constructor (rank 0)
//   PPole() : Matrix<double>::Matrix() {}

  /// Constructor 
  PPole(complex<double>& beta, vector<double>& g_i, double mass_0); 

  /// Copy Constructor
  // PPole(const PPole &theCopy);

  /// Destructor
  virtual ~PPole();

  // operators:


  // functions:

  virtual void evalMatrix(const double mass);

  void updateBeta (complex<double> beta) {_beta=beta;}

  
protected:
  complex<double> _beta;
};
//_____________________________________________________________________________
