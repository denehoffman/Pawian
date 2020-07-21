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

// KPoleBarrier class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KPoleBarrier.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include "PwaDynamics/KPole.hh"

#include <vector>
#include <memory>

class AbsPhaseSpace;

class KPoleBarrier : public KPole {

public:
  /// Constructor 
  KPoleBarrier(vector<double>& g_i, double mass_0, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs,
	       int orbMom, bool truncatedBarrier=false); 

  /// Destructor
  virtual ~KPoleBarrier();

  virtual void evalMatrix(const double mass, Spin OrbMom=0);
  virtual void updatePoleMass (double newPoleMass);

  
protected:
  vector<std::shared_ptr<AbsPhaseSpace> > _phpVecs;

};

