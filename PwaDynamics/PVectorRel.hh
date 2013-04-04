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
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>

class PPole;
class AbsPhaseSpace;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class PVectorRel : public KMatrixBase {

public:

  /// Constructor 
  PVectorRel(vector<boost::shared_ptr<PPole> > Ppoles, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs); 
  PVectorRel(vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs);

  /// Destructor
  virtual ~PVectorRel();

  virtual void evalMatrix(const double mass);

  virtual void updateBeta(int i, complex<double> beta);

protected:
  vector<boost::shared_ptr<PPole> > _Ppoles;
};
//_____________________________________________________________________________


