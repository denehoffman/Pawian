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

// FVectorPiPiS class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file FVectorPiPiS.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/PVectorRel.hh"
#include "PwaDynamics/PVectorSlowCorRel.hh"
#include <iostream>
#include <vector>
#include <memory>

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class FVectorPiPiS : public FVector {

public:

  /// Constructor 
  FVectorPiPiS(); 

  /// Destructor
  virtual ~FVectorPiPiS();
  virtual std::shared_ptr<PVectorSlowCorRel> pVectorSlowCorrRel() {return _pVectorCor;}
  virtual void updateFprod (int i, complex<double> fProd){_pVectorCor->updateFprod (i, fProd);}
  virtual void updateS0prod (double s0prod) {_pVectorCor->updateS0prod(s0prod);}
protected:
  std::shared_ptr<PVectorSlowCorRel> _pVectorCor;
};
//_____________________________________________________________________________


