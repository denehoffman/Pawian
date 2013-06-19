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

// Flatte class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file Flatte.h
//_____________________________________________________________________________

#include <iostream>
#include <vector>
#include <complex>
#include <utility>
#include <memory>

class Particle;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class Flatte {

public:

  /// Constructor 
  Flatte(std::pair<Particle*, Particle*>& decPair1, std::pair<Particle*, Particle*>& decPair2); //decPair1: decay particles
                                                                                                //decPair2: second decay channel    
  Flatte(std::pair <const double, const double>& massPair1, std::pair <const double, const double>& massPair2);
 
  /// Destructor
  ~Flatte();

  complex<double> calcFirstChannel(double currentMass, double mass0, double g1, double g2);
  complex<double> calcSecondChannel(double currentMass, double mass0, double g1, double g2);

protected:
  double _mass11;
  double _mass12;
  double _mass21;
  double _mass22; 
};



