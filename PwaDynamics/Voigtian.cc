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

#include "PwaDynamics/Voigtian.hh"
#include "Utils/Faddeeva.hh"

#include <complex>
using std::complex;

Voigtian::Voigtian() {
}

Voigtian::~Voigtian(){
}

double Voigtian::calc(double currentMass, double mass0, double gamma, double sigma){
  double result=0.;
  if ((sigma < 0. || gamma < 0.) || (fabs(sigma)<1e-20 && fabs(gamma)<1.e-20)) {
    return result;  // Not meant to be for those who want to be thinner than 0
  }
  if (fabs(sigma) < 1.e-20){
    sigma=1.e-12;
  }
  if (fabs(gamma) < 1.e-20){
    gamma=1.e-12;
  }
  
  double denom=sqrt(2.)*sigma;  
  double realZ=(currentMass-mass0)/denom;
  double imagZ=gamma/(2.*denom);
  complex<double> complZ(realZ,imagZ);  

  result=sqrt(2.*M_PI)/4.*gamma/sigma*Faddeeva::w(complZ).real();
  return result;
}

