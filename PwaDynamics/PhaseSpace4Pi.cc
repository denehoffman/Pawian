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

// this paramererizaton of this 4pi phasespce factor has been taken over from 
// the Laura++ software package:
// 
// Copyright University of Warwick 2008 - 2013.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// Authors:
// Thomas Latham
// John Back
// Paul Harrison
//
// Calculate the 4pi phase space factor. This uses a 6th-order polynomial 
// parameterisation that approximates the multi-body phase space double integral
// defined in Eq 4 of the A&S paper hep-ph/0204328. This form agrees with the 
// BaBar model (another 6th order polynomial from s^4 down to 1/s^2), but avoids the
// exponential increase at small values of s (~< 0.1) arising from 1/s and 1/s^2.
// Eq 4 is evaluated for each value of s by assuming incremental steps of 1e-3 GeV^2 
// for s1 and s2, the invariant energy squared of each of the di-pion states,
// with the integration limits of s1 = (2*mpi)^2 to (sqrt(s) - 2*mpi)^2 and
// s2 = (2*mpi)^2 to (sqrt(s) - sqrt(s1))^2. The mass M of the rho is taken to be
// 0.775 GeV and the energy-dependent width of the 4pi system 
// Gamma(s) = gamma_0*rho1^3(s), where rho1 = sqrt(1.0 - 4*mpiSq/s) and gamma_0 is 
// the "width" of the 4pi state at s = 1, which is taken to be 0.3 GeV 
// (~75% of the total width from PDG estimates of the f0(1370) -> 4pi state).
// The normalisation term rho_0 is found by ensuring that the phase space integral
// at s = 1 is equal to sqrt(1.0 - 16*mpiSq/s). Note that the exponent for this 
// factor in hep-ph/0204328 is wrong; it should be 0.5, i.e. sqrt, not n = 1 to 5.
// Plotting the value of this double integral as a function of s can then be fitted
// to a 6th-order polynomial (for s < 1), which is the result used below


#include "PwaDynamics/PhaseSpace4Pi.hh"
#include "Utils/PawianConstants.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

PhaseSpace4Pi::PhaseSpace4Pi():
  AbsPhaseSpace()
  ,_fourPiFactor1(16.0*PawianConstants::mPiSq)
{
}

PhaseSpace4Pi::~PhaseSpace4Pi(){

}

complex<double> PhaseSpace4Pi::factor(const double mass){
  double mass_sqr=mass*mass;
  complex<double> result(0.,0.);

  if (fabs(mass_sqr) < 1e-10) {return result;}

  if (mass_sqr <= 1.0) {
    double rhoTerm = ((1.07885*mass_sqr + 0.13655)*mass_sqr - 0.29744)*mass_sqr - 0.20840;
    rhoTerm = ((rhoTerm*mass_sqr + 0.13851)*mass_sqr - 0.01933)*mass_sqr + 0.00051;
    // For some values of s (below 2*mpi), this term is a very small 
    // negative number. Check for this and set the rho term to zero.
    if (rhoTerm < 0.0) {rhoTerm = 0.0;}
    result= complex<double>( rhoTerm, 0. );
  } else {
    result= complex<double>( sqrt(1.0 - (_fourPiFactor1/mass_sqr)), 0. );
  }  

  return result;
}

complex<double> PhaseSpace4Pi::factor(const complex<double> mass){
  complex<double> mass_sqr=mass*mass;
  complex<double> result(0.,0.);
  if( norm(mass_sqr) <= 1. ){
    complex<double> rhoTerm = ((1.07885*mass_sqr + 0.13655)*mass_sqr - 0.29744)*mass_sqr - 0.20840;
    rhoTerm = ((rhoTerm*mass_sqr + 0.13851)*mass_sqr - 0.01933)*mass_sqr + 0.00051;
    result = rhoTerm;
    }
  else result = sqrt(1.0 - (_fourPiFactor1/mass_sqr));
  CorrectForChosenSign(result, result);
  return result;
}

complex<double> PhaseSpace4Pi::breakUpMom(const double mass){
   return factor(mass)*mass/2.;
}

complex<double> PhaseSpace4Pi::breakUpMom(const complex<double> mass){
  complex<double> result = factor(mass)*mass/2.;
  return result;
   // Alert << "PhaseSpace4Pi does currently not support complex masses" << endmsg;
   // exit(EXIT_FAILURE);
   // return 0;
}

