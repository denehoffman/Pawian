//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universit??t Bochum 				  //
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

#include "PwaDynamics/PhaseSpaceIsobarReid.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"

PhaseSpaceIsobarReid::PhaseSpaceIsobarReid(double mass1, double mass2):
  PhaseSpaceIsobar(mass1, mass2)
{
}

PhaseSpaceIsobarReid::~PhaseSpaceIsobarReid(){

}

complex<double> PhaseSpaceIsobarReid::factor(const double mass){
   return PawianQFT::phaseSpaceFacReid(mass, _mass1, _mass2);
}

complex<double> PhaseSpaceIsobarReid::breakUpMom(const double mass){
  return breakupMomQ(mass,_mass1, _mass2);
}

complex<double> PhaseSpaceIsobarReid::factor(const complex<double> mass){
  complex<double> rho = PawianQFT::phaseSpaceFacReid(mass, _mass1, _mass2);
  complex<double> mom=rho*mass/2.;
  //  CorrectForChosenSign(mom, rho);
  CorrectForChosenSign(rho, rho);
  return rho;
}

complex<double> PhaseSpaceIsobarReid::breakUpMom(const complex<double> mass){
  complex<double> q = breakupMomQ(mass,_mass1, _mass2);
  CorrectForChosenSign(q, q);
  return q;
}

void PhaseSpaceIsobarReid::CorrectForChosenSign(complex<double>& breakUpMom, complex<double>& toChange){
    if( (_bumImPartSign > 0 && breakUpMom.real()>0.) || (_bumImPartSign < 0 && fabs(breakUpMom.real())<1.e-10 ) ){
      toChange *= -1; // for Reid php factors
    }
  }

