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

// KMatrixSlowAdlerCor class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KMatrixSlowAdlerCorRel.h
//_____________________________________________________________________________

#include "PwaDynamics/KMatrixBase.hh"
#include <vector>
#include <memory>
#include <boost/multi_array.hpp>

class KPole;
class AbsPhaseSpace;

typedef boost::multi_array< double, 2> array_type_2d;

class KMatrixSlowAdlerCorRel : public KMatrixBase {

public:

  /// Constructor 
  KMatrixSlowAdlerCorRel(vector<std::shared_ptr<KPole> > Kpoles,
			 vector<std::shared_ptr<AbsPhaseSpace> > phpVecs,
			 std::shared_ptr<array_type_2d> fscatPtr, double s0_scat,
			 double sAdler0=-0.15, double sAdler=1.); 
  KMatrixSlowAdlerCorRel(int dim);
  /// Destructor
  virtual ~KMatrixSlowAdlerCorRel();

  virtual void evalMatrix(const double mass, Spin OrbMom=0);
  virtual void evalMatrix(const complex<double> mass, Spin OrbMom=0);
  virtual std::shared_ptr<array_type_2d> fScatProd() {return _fScatPtr;}
  virtual double s0Scat() {return _s0Scat;}

protected:
  std::shared_ptr<array_type_2d> _fScatPtr;
  double _s0Scat; 
  double _sAdler0;
  double _sAdler;

private:
  template<typename MassType>
  void evalMatrixTemplate(const MassType mass, Spin OrbMom=0);
};



