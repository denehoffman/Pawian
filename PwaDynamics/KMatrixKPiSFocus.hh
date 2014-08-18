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

//KMatrixKPiSFocus class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @fileKMatrixKPiSFocus.h
//_____________________________________________________________________________

#include "PwaDynamics/KMatrixBase.hh"
#include <iostream>
#include <vector>
#include <memory>
#include <boost/multi_array.hpp>

class KPole;
class AbsPhaseSpace;

typedef boost::multi_array< double, 2> array_type_2d;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class KMatrixKPiSFocus : public KMatrixBase {

public:

  /// Constructor 
  KMatrixKPiSFocus(int isospin2=1); 
  /// Destructor
  virtual ~KMatrixKPiSFocus();

  virtual void evalMatrix(const double mass);
  virtual void evalMatrix(const complex<double> mass);
  virtual double sNorm() {return _sNorm;}

protected:
  std::shared_ptr<array_type_2d> _aScatPtr;
  std::shared_ptr<array_type_2d> _bScatPtr;
  std::shared_ptr<array_type_2d> _cScatPtr;
  double _sAdler0;
  const double _sNorm;

  void init2IsoSpin1();
  void init2IsoSpin3();

private:
  template<typename MassType>
  void evalMatrixTemplate(const MassType mass);
};
//_____________________________________________________________________________


