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

#include <vector>
#include <complex>
#include <memory>

class AbsPhaseSpace;

class KPole : public Matrix< complex<double> > {

public:

  /// Constructor 
  KPole(vector<double>& g_i, double mass_0);
  KPole(vector<double>& g_i, double mass_0, int numRow, int numCol); 

  /// Destructor
  virtual ~KPole();

  virtual void evalMatrix(const double mass, Spin OrbMom=0);
  virtual void evalMatrix(const complex<double> mass, Spin OrbMom=0);

  virtual double poleMass() {return _poleMass;}
  virtual vector<double> gFactors() {return _g_i;}

  virtual void updatePoleMass (double newPoleMass);
  virtual void updategFactors (vector<double>& newg_i) {_g_i=newg_i;}

  virtual std::vector< complex<double> > barrierFactors() {return _barrierFactor;}
  
protected:
  vector<double> _g_i;
  double _poleMass;
  complex<double> _poleMassCompl;

  int _orbMom;
  std::vector< complex<double> > _breakUpM0;
  std::vector< complex<double> > _barrierFactor;
  bool _truncatedBarrier;

private:

};

