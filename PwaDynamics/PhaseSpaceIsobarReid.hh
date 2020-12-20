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

// PhaseSpaceIsobarReid class definition file. -*- C++ -*-
// Copyright 2016 Bertram Kopf
// phase space factor which fulfils the analytic continuation and unitarity

#pragma once 
//_____________________________________________________________________________
// @file PhaseSpaceIsobarReid.h
//_____________________________________________________________________________

#include "PwaDynamics/PhaseSpaceIsobar.hh"

#include <complex>
using std::complex;

class PhaseSpaceIsobarReid: public PhaseSpaceIsobar {

public:

  /// Constructor 
  PhaseSpaceIsobarReid(double mass1, double mass2); 

  /// Destructor
  virtual ~PhaseSpaceIsobarReid();

  virtual complex<double> factor(const double mass);
  virtual complex<double> breakUpMom(const double mass);
  virtual complex<double> factor(const complex<double> mass);
  virtual complex<double> breakUpMom(const complex<double> mass);
  virtual complex<double> ChewM(const double mass, int orbMom=0);
  virtual complex<double> ChewM(const complex<double> mass, int orbMom=0);

protected:
  void CorrectCMForChosenSign(complex<double>& breakUpMom, complex<double>& toChange);

};



