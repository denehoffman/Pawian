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

//KMatrixPiPiS class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KMatrixPiPiS.h
//_____________________________________________________________________________

#include "PwaDynamics/KMatrixSlowAdlerCorRel.hh"
#include <iostream>
#include <vector>
#include <memory>
#include <boost/multi_array.hpp>


typedef boost::multi_array< double, 2> array_type_2d;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class KMatrixPiPiS : public KMatrixSlowAdlerCorRel {

public:

  /// Constructor 
 KMatrixPiPiS(); 

  /// Destructor
  virtual ~KMatrixPiPiS();

  // virtual void evalMatrix(const double mass);

protected:
  void initASParam1900();
  double _oldMass;
private:
};
//_____________________________________________________________________________


