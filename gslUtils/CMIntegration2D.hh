//************************************************************************//
//									  //
//  Copyright 2023 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Meike Kuessner (mkuessner@ep1.rub.de)		  //
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

// Voitian class definition file. -*- C++ -*-
// Copyright 2023 Bertram Kopf

#pragma once 

#include <complex>
#include <stddef.h>
#include <gsl/gsl_integration.h>

#include "gslUtils/CMIntegration.hh"


//_____________________________________________________________________________
// @file CMIntegration2D.hh
//_____________________________________________________________________________

class CMIntegration2D :  public CMIntegration{

public:
  /// Constructor 
  CMIntegration2D(double mpole1, double fpole1, double m11, double m21,double mpole2, double fpole2, double m12, double m22); 
 
  /// Destructor
  ~CMIntegration2D();

  virtual void doFit2(gsl_function& F, double& result, double& resulterr, std::string fitName);
  void setCMparams2(CMunstable_params& theParams);
  
protected:
  static CMunstable_params _CMunstable_params2;
  static std::complex<double> _currentS2;
  static double _currentsprimeprime;
  static std::complex<double> _currentOuterResult;
  static std::complex<double> _currentInnerResult;
  static bool _calcRealPartInner;
  double _integLowerBorder2;
  
private:  
};



