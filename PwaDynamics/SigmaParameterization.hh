//************************************************************************//
//									  //
//  Copyright 2024 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// SigmaParameterization class definition file. -*- C++ -*-
// Copyright 2024 Bertram Kopf

// Fixed parametrization for the sigma pole based on
// Phys. Lett B 598 (204) 149-158

#pragma once 
#include <complex>
using std::complex;

//_____________________________________________________________________________
// @file SigmaParameterization.h
//_____________________________________________________________________________

class SigmaParameterization {

public:
  /// Constructor
  SigmaParameterization();  
  SigmaParameterization(double mPole, double g4pi, double b1, double b2, double a); 
 
  /// Destructor
  ~SigmaParameterization();

  complex<double> calc(double currentMass);
  complex<double> calc(double currentMass, std::complex<double> gSigma);
  complex<double> calcT(double currentMass);

  const double barePoleMass() const {return _mPole;}
  void setBarePoleMass(double barePoleMass);
  const double b1Param() const {return _b1;}
  void setb1Param(double b1) {_b1=b1;}
  const double b2Param() const {return _b2;}
  void setb2Param(double b2) {_b2=b2;}
  const double aParam() const {return _a;}
  void setaParam(double a) {_a=a;}
  const double g4piParam() const {return _g4pi;}
  void setg4piParam(double g4pi) {_g4pi=g4pi;}
  
protected:
  std::complex<double> php4pi(double currentMass2);
  double _mPole;
  double _g4pi;
  double _b1;
  double _b2;
  double _a;
  const double _sAdler;
  double _mPoleSqr;
  double _g2piDenom;
  std::complex<double> _rho2pi_mPole;
  std::complex<double> _rho4pi_mPole;
};



