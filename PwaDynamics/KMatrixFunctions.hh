//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// Copyright 2017 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KMatrixFunctions.hh
//_____________________________________________________________________________

#include <iostream>
#include <complex>
#include <utility>
#include <memory>
#include <map>

#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaDynamics/BarrierFactor.hh"


class TMatrixRel;
//_____________________________________________________________________________
//_____________________________________________________________________________

namespace KMatrixFunctions { 
  double twoDeltaArgand(std::shared_ptr<TMatrixRel> tMatr, unsigned int projectionIndex, double mass);
  double deltaArgand(std::shared_ptr<TMatrixRel> tMatr, unsigned int projectionIndex, double mass);
  unsigned int noOfPhaseRotationsArgand(std::shared_ptr<TMatrixRel> tMatr, unsigned int projectionIndex, double mass);
  void fillRotationArgandMap(std::shared_ptr<TMatrixRel> tMatr, unsigned int projectionIndex, std::map<unsigned int, double>& toFill, double massMax=4.0);
}; // namespace KMatrixFunctions



