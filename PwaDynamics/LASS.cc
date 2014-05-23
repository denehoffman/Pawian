//************************************************************************//
//									  //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#include "qft++/relativistic-quantum-mechanics/Utils.hh"

#include "PwaDynamics/LASS.hh"

complex<double>  LASS::K0star_1430(double currentMass, double m0, double gammaM, double aLASS, double rLASS, double BLASS, double phiB, double RLASS, double phiR){
  complex<double> result(0.,0.);
  if(currentMass>1.6) return result;

  double KMass=0.49367;
  double piMass=0.13498;

  double breakupQ=breakupMomQ(currentMass, KMass, piMass).real();

  //calculate the background phase motion
  double cot_deltaB = 1.0/(aLASS*breakupQ) + 0.5*rLASS*breakupQ;
  double deltaB = atan( 1.0/cot_deltaB);
  double totalB = deltaB + phiB ;

  //calculate the resonant phase motion
  double deltaR = atan((m0*gammaM/(m0*m0 - currentMass*currentMass)));
  double totalR = deltaR + phiR ;

  //sum them up
  complex<double> bkgB = complex<double>(BLASS*sin(totalB),0)*complex<double>(cos(totalB),sin(totalB));
  complex<double> resT = complex<double>(rLASS*sin(deltaR),0)*complex<double>(cos(totalR),sin(totalR))*complex<double>(cos(2*totalB),sin(2*totalB));
  result = bkgB + resT;      

  return result;
}

