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
//									  //
//************************************************************************//

// KMatrixOmnes class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KMatrixOmnes.h
//_____________________________________________________________________________

#include "PwaDynamics/KMatrixNonRel.hh"
#include <vector>
#include <memory>

class KPoleOmnes;
class AbsPhaseSpace;

class KMatrixOmnes : public KMatrixNonRel {

public:

  /// Constructor 
  KMatrixOmnes(vector<std::shared_ptr<KPole> > KpolesOmnes, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs); 

  /// Destructor
  virtual ~KMatrixOmnes();

  virtual void evalMatrix(const double mass, Spin OrbMom=0);
  virtual void evalMatrix(const complex<double> mass, Spin OrbMom=0);
protected:

private:
  template<typename MassType>
  void evalMatrixTemplate(const MassType mass, Spin OrbMom=0);
};



