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

// PhaseSpaceIsobar class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 
//_____________________________________________________________________________
// @file PhaseSpaceIsobar.h
//_____________________________________________________________________________

#include <iostream>
#include <vector>
#include "PwaDynamics/AbsPhaseSpace.hh"

using namespace std;

//_____________________________________________________________________________
//_____________________________________________________________________________

class PhaseSpaceIsobar: public AbsPhaseSpace {

public:

  /// Constructor 
  PhaseSpaceIsobar(double mass1, double mass2); 

  /// Destructor
  virtual ~PhaseSpaceIsobar();

  // operators:


  // functions:

  virtual complex<double> factor(const double mass);
  virtual complex<double> breakUpMom(const double mass);
  virtual complex<double> factor(const complex<double> mass);
  virtual complex<double> breakUpMom(const complex<double> mass);

protected:

private:
  double _mass1;
  double _mass2;
};
//_____________________________________________________________________________


