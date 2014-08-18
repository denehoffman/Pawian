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

#include "PwaDynamics/PhaseSpace4Pi.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

PhaseSpace4Pi::PhaseSpace4Pi():
  AbsPhaseSpace()
  , _piMass(0.13957)
{
}

PhaseSpace4Pi::~PhaseSpace4Pi(){

}

complex<double> PhaseSpace4Pi::factor(const double mass){
  double mass_sqr=mass*mass;
  complex<double> result(0.,0.);
  if( mass_sqr <= 1 ){
      double real = 1.2274 + .00370909 / ( mass_sqr * mass_sqr ) - .111203 / mass_sqr - 6.39017 * mass_sqr + 
	 16.8358*mass_sqr*mass_sqr - 21.8845*mass_sqr*mass_sqr*mass_sqr + 11.3153*mass_sqr*mass_sqr*mass_sqr*mass_sqr;

      double cont32 = sqrt(1.0-(16.0*_piMass*_piMass));
      result = complex<double>( cont32 * real, 0 );
    }
  else result = complex<double>( sqrt( 1. - 16. * _piMass * _piMass / mass_sqr ), 0. );
  return result;
}

complex<double> PhaseSpace4Pi::factor(const complex<double> mass){
   Alert << "PhaseSpace4Pi does currently not support complex masses" << endmsg;
   exit(EXIT_FAILURE);
   return 0;
}

complex<double> PhaseSpace4Pi::breakUpMom(const double mass){
   return factor(mass)*mass/2.;
}

complex<double> PhaseSpace4Pi::breakUpMom(const complex<double> mass){
   Alert << "PhaseSpace4Pi does currently not support complex masses" << endmsg;
   exit(EXIT_FAILURE);
   return 0;
}

