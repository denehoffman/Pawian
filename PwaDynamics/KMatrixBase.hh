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

// KMatrixBase class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KMatrixBase.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include <iostream>
#include <vector>
#include <memory>
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class KPole;


using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class KMatrixBase : public Matrix< complex<double> > {

public:

  /// Constructor 
  KMatrixBase(vector<std::shared_ptr<KPole> > Kpoles, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs); 
  KMatrixBase(vector<std::shared_ptr<AbsPhaseSpace> > phpVecs, int numCols, int numRows);
  KMatrixBase(int numCols, int numRows); 
  /// Destructor
  virtual ~KMatrixBase();

  virtual void evalMatrix(const double mass, Spin OrbMom=0);
  virtual void evalMatrix(const complex<double> mass, Spin OrbMom=0) {return;}
  virtual vector<std::shared_ptr<AbsPhaseSpace> > phaseSpaceVec() {return _phpVecs;}
  virtual vector<std::shared_ptr<KPole> > kpoles() {return _KPoles;}
  virtual void updateBgTerms(unsigned int order, unsigned int row,  unsigned int column, double theVal);
  virtual void updates0Adler(double s0Adler) {_s0Adler=s0Adler;}
  virtual void updatesnormAdler(double snormAdler) {_snormAdler=snormAdler;}

  void SetBumImPartSigns(std::vector<double> signs);

protected:
  vector<std::shared_ptr<KPole> > _KPoles;
  vector<std::shared_ptr<AbsPhaseSpace> > _phpVecs;
  vector< vector< vector<double> > > _bgTerms;
  unsigned int _orderBg;
  double _s0Adler;
  double _snormAdler;
};
//_____________________________________________________________________________


