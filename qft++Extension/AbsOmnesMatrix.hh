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

// AbsOmnesMatrix class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#pragma once 
//_____________________________________________________________________________
// @file PhaseSpaceIsobar.h
//_____________________________________________________________________________
#include <vector>
#include <memory>
#include <complex>
#include <boost/multi_array.hpp>

#include "qft++Extension/AbsComplexFunction.hh"



using std::complex;

class AbsOmnesMatrix{

public:

  // Constructor 
  AbsOmnesMatrix(){;} 

  /// Destructor
  virtual ~AbsOmnesMatrix() {;}
  virtual const std::string name()=0;
  virtual boost::multi_array< std::shared_ptr<AbsComplexFunction> , 2> omnesMatrix()=0;
  virtual boost::multi_array< std::shared_ptr<AbsComplexFunction> , 2> selfEnergyMatrix()=0;
  
protected:

};


