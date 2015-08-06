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

// PVectorSlowSlowCorRel class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file PVectorSlowCorRel.h
//_____________________________________________________________________________

#include "PwaDynamics/PVectorRel.hh"
#include <iostream>
#include <vector>
#include <memory>
//#include <boost/multi_array.hpp>

class PPole;
class AbsPhaseSpace;

//typedef boost::multi_array< complex<double>, 2> array_type_2dc;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class PVectorSlowCorRel : public PVectorRel {

public:

  /// Constructor 
  PVectorSlowCorRel(vector<std::shared_ptr<PPole> > Ppoles, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs, std::vector< complex<double> >& fProdVec, double s0prod);

  /// Destructor
  virtual ~PVectorSlowCorRel();

  virtual void evalMatrix(const double mass, Spin OrbMom=0);
  virtual void updateFprod (int i, complex<double> fProd);
  virtual void updateS0prod (double s0prod) {_s0prod=s0prod;}

protected:
  std::vector< complex<double> > _fProdVec;
  double _s0prod; 
};
//_____________________________________________________________________________


