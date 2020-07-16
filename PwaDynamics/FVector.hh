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
#include <vector>
#include <memory>

#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaDynamics/PVectorRel.hh"

class AbsPhaseSpace;
class KMatrixBase;

class FVector : public Matrix< complex<double> > {

public:

  FVector(std::shared_ptr<KMatrixBase> Kmatrix, std::shared_ptr<PVectorRel> Pvector);
  FVector(int numRows); 
  virtual ~FVector();

  virtual void evalMatrix(const double mass, Spin OrbMom=0);
  virtual void evalMatrix(const complex<double> mass, Spin OrbMom=0);
  virtual complex<double> evalProjMatrix(const double mass, int index, Spin OrbMom=0);
  virtual std::shared_ptr<KMatrixBase> kMatrix(){return _Kmatrix;}
  virtual std::shared_ptr<KMatrixBase> pVector(){return _Pvector;}
  void SetBumImPartSigns(std::vector<double> signs);

protected:
  std::shared_ptr<KMatrixBase> _Kmatrix; 
  std::shared_ptr<PVectorRel> _Pvector;
  complex<double> _imagCompl;
  IdentityMatrix< complex<double> > _idMatrix;
  Matrix< complex<double> > _CMMatrix;
  vector<std::shared_ptr<AbsPhaseSpace> > _phpVec;
  double _cSign;
};

