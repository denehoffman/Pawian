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

// KPole class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KPole.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <iostream>
#include <vector>
#include <complex>
#include <memory>
using namespace std;

class AbsPhaseSpace;
//_____________________________________________________________________________
//_____________________________________________________________________________

class KPole : public Matrix< complex<double> > {

public:

  // create/copy/destroy:

  /// Default Constructor (rank 0)
//   KPole() : Matrix<double>::Matrix() {}

  /// Constructor 
  KPole(vector<double>& g_i, double mass_0);
  KPole(vector<double>& g_i, double mass_0, int numRow, int numCol); 

  /// Copy Constructor
  // KPole(const KPole &theCopy);

  /// Destructor
  virtual ~KPole();

  // operators:


  // functions:

  virtual void evalMatrix(const double mass, Spin OrbMom=0);
  virtual void evalMatrix(const complex<double> mass, Spin OrbMom=0);

  virtual double poleMass() {return _poleMass;}
  virtual vector<double> gFactors() {return _g_i;}

  //  virtual void updatePoleMass (double newPoleMass) {_poleMass=newPoleMass;}
  virtual void updatePoleMass (double newPoleMass);
  virtual void updategFactors (vector<double>& newg_i) {_g_i=newg_i;}
  
protected:
  vector<double> _g_i;
  double _poleMass;
  complex<double> _poleMassCompl;

private:

};
//_____________________________________________________________________________
