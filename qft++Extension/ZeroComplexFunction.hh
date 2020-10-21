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
//	      							  //
//************************************************************************//

// ZeroComplexFunction class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#pragma once 
//_____________________________________________________________________________
// @file PhaseSpaceIsobar.h
//_____________________________________________________________________________
#include <vector>
#include <memory>
#include <complex>

#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++Extension/PawianUtils.hh"



using std::complex;

class ZeroComplexFunction : public AbsComplexFunction{

public:

  /// Constructor 
  ZeroComplexFunction(){;} 

  /// Destructor
  virtual ~ZeroComplexFunction(){;}
  virtual complex<double> eval(double s) {return complex<double>(1.,0.);}
  virtual complex<double> eval(complex<double> s) {return complex<double>(1.,0.);}
  
protected:
  
};


