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

// PVectorRel class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file PVectorRel.h
//_____________________________________________________________________________

#include "PwaDynamics/KMatrixBase.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include <iostream>
#include <vector>
#include <memory>

class PPole;
class AbsPhaseSpace;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class PVectorRel : public KMatrixBase {

public:

  /// Constructor 
  PVectorRel(vector<std::shared_ptr<PPole> > Ppoles, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs); 
  PVectorRel(vector<std::shared_ptr<AbsPhaseSpace> > phpVecs);

  /// Destructor
  virtual ~PVectorRel();

  virtual void evalMatrix(const double mass, Spin OrbMom=0);

  virtual void updateBeta(int i, complex<double> beta);
  virtual void updatePoleMass(int i, double mass);
  virtual void updategFactors(int i, vector<double>& newg_i);
  virtual void updateBgTerms(unsigned int order, unsigned int channel, double theVal) {return;}
  virtual void printElements();

protected:
  vector<std::shared_ptr<PPole> > _Ppoles;
  unsigned int _nOfChannels;
};
//_____________________________________________________________________________


