//************************************************************************//
//									  //
//  Copyright 2023 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Meike Kuessner (mkuessner@ep1.rub.de)		  //
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

// Voitian class definition file. -*- C++ -*-
// Copyright 2023 Bertram Kopf

#pragma once 

#include <complex>
#include <stddef.h>

#include "gslUtils/CMIntegration.hh"


//_____________________________________________________________________________
// @file CMIntegrationReid.hh
//_____________________________________________________________________________

class CMIntegrationReid :  public CMIntegration {

public:
  /// Constructor 
  CMIntegrationReid(double mpole, double fpole, double mu, double m1, double m2); 
 
  /// Destructor
  ~CMIntegrationReid();

  virtual void integrate(std::complex<double> s, std::complex<double>& result, std::complex<double>& resulterr);

protected:

private:
  static double FIntWrapper(double x, void * params);
  static double FIntWrapperNormCondition(double x, void * params);
  static std::complex<double> Ctilde(std::complex<double> s, double sprime);
  static std::complex<double> Sigma(double sprime, double m1, double m2);
  static double dsNorm(double sprime, double m1, double m2, double mpole, double fpole);
  static double dsImag(double sprime, double m1, double m2, double mpole, double fpole);
};



