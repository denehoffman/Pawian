 //************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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
//	      							  //
//************************************************************************//

// PhaseSpaceList class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#pragma once 
//_____________________________________________________________________________
// @file PhaseSpaceIsobar.h
//_____________________________________________________________________________
#include <vector>
#include <memory>
#include <complex>

#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/AbsComplexFunction.hh"
#include "qft++Extension/PawianUtils.hh"



using std::complex;

class PhaseSpaceList{

public:

  /// Constructor 
  PhaseSpaceList(vector<std::shared_ptr<AbsPhaseSpace> > phpVecs); 

  /// Destructor
  virtual ~PhaseSpaceList();
  vector< std::shared_ptr<AbsPhaseSpace> > phpVecs(){return _phpVecs;}
  Matrix< std::shared_ptr<AbsComplexFunction> > omnesMatrix(){return _omnesMatr;}
  Matrix< std::shared_ptr<AbsComplexFunction> > selfEnergyMatrix(){return _selfEnergyMatr;}
  
protected:
  vector< std::shared_ptr<AbsPhaseSpace> > _phpVecs;
  Matrix< std::shared_ptr<AbsComplexFunction> > _omnesMatr;
  Matrix< std::shared_ptr<AbsComplexFunction> > _selfEnergyMatr;
};


