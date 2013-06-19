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

// FVector class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file FVector.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include "PwaDynamics/KMatrixBase.hh"
#include "PwaDynamics/PVectorRel.hh"
#include <iostream>
#include <vector>
#include <memory>

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class FVector : public Matrix< complex<double> > {

public:

  /// Constructor 
  FVector(std::shared_ptr<KMatrixBase> Kmatrix, std::shared_ptr<PVectorRel> Pvector);
  FVector(int numRows); 

  /// Destructor
  virtual ~FVector();

  virtual void evalMatrix(const double mass);
  virtual void updateBeta(int i, complex<double> beta) {_Pvector->updateBeta(i, beta);}
  //  virtual void updateFprod (int i, complex<double> fProd) {;}
  //  virtual void updateS0prod (double s0prod) {;}
  virtual std::shared_ptr<KMatrixBase> kMatrix(){return _Kmatrix;}
protected:
  std::shared_ptr<KMatrixBase> _Kmatrix; 
  std::shared_ptr<PVectorRel> _Pvector; 
};
//_____________________________________________________________________________


