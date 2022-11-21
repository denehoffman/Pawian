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

struct CMunstable_params {
  double _mPole;
  double _fPole;
  double _mu;
  double _m1;
  double _m2;

  CMunstable_params(double mpole, double fpole, double mu, double m1, double m2) :
    _mPole(mpole)
    ,_fPole(fpole)
    ,_mu(mu)
    ,_m1(m1)
    ,_m2(m2){
  }
};

//_____________________________________________________________________________
// @file CMIntegration.hh
//_____________________________________________________________________________

class CMIntegration {

public:
  /// Constructor 
  CMIntegration(double mpole, double fpole, double mu, double m1, double m2); 
 
  /// Destructor
  ~CMIntegration();

  std::complex<double> integrate(std::complex<double> s);

protected:

private:
  static double FIntWrapper(double x, void * params);
  static std::complex<double> Ctilde(std::complex<double> s, double sprime);
  static std::complex<double> Sigma(double sprime, double m1, double m2);
  static double dsNorm(std::complex<double> s, double m1, double m2, double mpole, double fpole);

  static CMunstable_params _CMunstable_params;
  static std::complex<double> _currentS;
  static bool _calcRealPart;
};



