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

// Copyright 2014 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file RadMultipoleFormFactor.hh
//_____________________________________________________________________________

#include <iostream>
#include <complex>
#include <utility>
#include <memory>

#include "qft++/relativistic-quantum-mechanics/Utils.hh"

//_____________________________________________________________________________
//_____________________________________________________________________________

namespace RadMultipoleFormFactor { 
  complex<double> PureM1(double massA, double massB, double currentMassB, double Egamma);

}; // namespace RadMultipoleFormFactor
