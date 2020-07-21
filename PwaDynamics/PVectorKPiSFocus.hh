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

// PVectorKPiSFocus class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file PVectorKPiSFocus.h
//_____________________________________________________________________________

#include "PwaDynamics/PVectorRel.hh"
#include "PwaDynamics/KMatrixKPiSFocus.hh"

#include <vector>
#include <memory>

class PVectorKPiSFocus : public PVectorRel {

public:
  /// Constructor 
  PVectorKPiSFocus(std::shared_ptr<KMatrixKPiSFocus> kMatrix);

  /// Destructor
  virtual ~PVectorKPiSFocus();

  virtual void evalMatrix(const double mass, Spin OrbMom=0);
  virtual void updateAprod (int i, double aProd);
  virtual void updateBprod (int i, double bProd);
  virtual void updateCprod (int i, double cProd);
  virtual void updatePhaseprod (int i, double phaseProd);

protected:
  std::vector< double > _aProdVec;
  std::vector< double > _bProdVec;
  std::vector< double > _cProdVec;
  std::vector< double > _phaseProdVec;
  std::shared_ptr<KMatrixKPiSFocus> _kMatrix;
  void init2IsoSpin1();
  void init2IsoSpin3();
};


