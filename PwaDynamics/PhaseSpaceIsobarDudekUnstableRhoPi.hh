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

// PhaseSpaceIsobarDudekUnstableRhoPi class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf
// phase space factor which fulfils the analytic continuation and unitarity

#pragma once 
//_____________________________________________________________________________
// @file PhaseSpaceIsobarDudekUnstableRhoPi.h
//_____________________________________________________________________________

#include "PwaDynamics/PhaseSpaceIsobar.hh"

#include <complex>

class PhaseSpaceIsobarDudekUnstableRhoPi : public PhaseSpaceIsobar {

public:

  /// Constructor 
  PhaseSpaceIsobarDudekUnstableRhoPi(double mass1, double mass2); 

  /// Destructor
  virtual ~PhaseSpaceIsobarDudekUnstableRhoPi();

  virtual std::complex<double> factor(const double mass);
  virtual std::complex<double> breakUpMom(const double mass);
  virtual std::complex<double> factor(const std::complex<double> mass);
  virtual std::complex<double> breakUpMom(const std::complex<double> mass);
  virtual std::complex<double> ChewM(const double mass);
  virtual std::complex<double> ChewM(const std::complex<double> mass);

protected:
  void CorrectCMForChosenSign(std::complex<double>& breakUpMom, std::complex<double>& toChange);

};


