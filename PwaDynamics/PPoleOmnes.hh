//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// PPoleOmnes class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file PPoleOmnes.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaDynamics/PPole.hh"
//#include "PwaDynamics/AbsPhaseSpace.hh"

#include <vector>
#include <complex>
using std::complex;
#include <memory>

class AbsPhaseSpace;

class PPoleOmnes : public PPole {

public:
  /// Constructor 
  PPoleOmnes(complex<double>& beta, vector<double>& g_i, double mass_0); 

  /// Destructor
  virtual ~PPoleOmnes();

  virtual void evalMatrix(const double mass, Spin OrbMom=0);
  virtual void evalMatrix(const complex<double> mass, Spin OrbMom=0);
  
protected:
};
