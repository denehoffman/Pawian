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

Voigtian::Voigtian()
 {
 }

Voigtian::~Voigtian(){
}

complex<double> Voigtian::calc(double currentMass, double mass0, double gamma, double sigma){
  
  double denom=sqrt(2.)*sigma;  
  double realZ=(currentMass-mass0)/denom;
  double imagZ=gamma/(2.*denom);
  complex<double> complZ(realZ,imagZ);  

  complex<double> result(sqrt(2.*M_PI)/4.*gamma/sigma*Faddeeva::erf(complZ).real(), 0.);
  return result;
}

