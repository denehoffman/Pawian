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

// PhaseSpace4Pi class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 
//_____________________________________________________________________________
// @file PhaseSpace4Pi.h
//_____________________________________________________________________________

#include "PwaDynamics/AbsPhaseSpace.hh"
#include "Utils/PawianConstants.hh"
#include "qft++Extension/PawianUtils.hh"

#include <complex>
using std::complex;

class PhaseSpace4Pi: public AbsPhaseSpace {

public:

  /// Constructor 
  PhaseSpace4Pi(); 

  /// Destructor
  virtual ~PhaseSpace4Pi();

  virtual complex<double> factor(const double mass, int orbMom=0);
  virtual complex<double> breakUpMom(const double mass);
  virtual complex<double> factor(const complex<double> mass, int orbMom=0);
  virtual complex<double> breakUpMom(const complex<double> mass);
  virtual complex<double> breakUpMomDefaultAS(const double mass) {
    return PawianQFT::breakupMomQDefaultAS(mass, 2.*PawianConstants::mPi, 2.*2.*PawianConstants::mPi);
  }
  virtual complex<double> breakUpMomDefaultAS(const complex<double> mass) {
    return PawianQFT::breakupMomQDefaultAS(mass, 2.*PawianConstants::mPi, 2.*2.*PawianConstants::mPi);
  }
  virtual double thresholdMass();

private:
  const double _fourPiFactor1;
};



