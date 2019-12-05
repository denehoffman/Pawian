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

// AbsPhaseSpace class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 
//_____________________________________________________________________________
// @file KMatrixBase.h
//_____________________________________________________________________________

#include <iostream>
#include <vector>
#include <complex>
#include <map>

#include "Utils/PawianConstants.hh"

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class AbsPhaseSpace {

public:

  /// Constructor 
  AbsPhaseSpace() : _bumImPartSign(1){;} 

  /// Destructor
  virtual ~AbsPhaseSpace(){;}

  // operators:


  // functions:

  virtual complex<double> factor(const double mass)=0;
  virtual complex<double> factor(const complex<double> mass)=0;
  virtual complex<double> breakUpMom(const double mass)=0;
  virtual complex<double> breakUpMom(const complex<double> mass)=0;
  virtual complex<double> breakUpMomDefaultAS(const double mass)=0;
  virtual complex<double> breakUpMomDefaultAS(const complex<double> mass)=0; 
  virtual complex<double> ChewM(const double mass) {return PawianConstants::i*factor(mass);}
  virtual complex<double> ChewM(const complex<double> mass) {return PawianConstants::i*factor(mass);}
  virtual double thresholdMass()=0;
  virtual void cacheFactors(const double mass) {return;}

  void SetBumImPartSign(double sign){_bumImPartSign = sign;}
  const std::string name() const {return _name;} 

protected:
  virtual void CorrectForChosenSign(complex<double>& breakUpMom, complex<double>& toChange){
   if((_bumImPartSign > 0 && breakUpMom.imag() < 0) ||
       (_bumImPartSign < 0 && breakUpMom.imag() > 0)){
      toChange *= -1;
    }
   }
  double _bumImPartSign;
  std::map<int, complex<double> > _CMCache; //resolution 100 keV

  std::string _name;

private:
  //  double _bumImPartSign;
};
//_____________________________________________________________________________


