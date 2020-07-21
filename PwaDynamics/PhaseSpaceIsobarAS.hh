//************************************************************************//
//									  //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// PhaseSpaceIsobarAS class definition file. -*- C++ -*-
// Copyright 2014 Julian Pychy

#pragma once 
//_____________________________________________________________________________
// @file PhaseSpaceIsobarAS.h
//_____________________________________________________________________________

#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++Extension/PawianUtils.hh"

#include <complex>

class PhaseSpaceIsobarAS: public AbsPhaseSpace {

public:

  /// Constructor 
  PhaseSpaceIsobarAS(double mass1, double mass2); 

  /// Destructor
  virtual ~PhaseSpaceIsobarAS();

  virtual std::complex<double> factor(const double mass);
  virtual std::complex<double> breakUpMom(const double mass);
  virtual std::complex<double> factor(const std::complex<double> mass);
  virtual std::complex<double> breakUpMom(const std::complex<double> mass);
  virtual std::complex<double> breakUpMomDefaultAS(const double mass) {
    return PawianQFT::breakupMomQDefaultAS(mass, _mass1, _mass2);
  }
  virtual std::complex<double> breakUpMomDefaultAS(const std::complex<double> mass) {
    return PawianQFT::breakupMomQDefaultAS(mass, _mass1, _mass2);
  }
  virtual double thresholdMass();

protected:
  double _mass1;
  double _mass2;

};


