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

#include "ErrLogger/ErrLogger.hh"
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
  //return PawianQFT::breakupMomQReid(mass,_mass1, _mass2);
  //  return PawianQFT::breakupMomQDefault(mass,_mass1, _mass2);
  return PawianQFT::breakupMomQDefaultAS(mass,_mass1, _mass2);
}

complex<double> PhaseSpaceIsobarReid::factor(const complex<double> mass){
  complex<double> rho = PawianQFT::phaseSpaceFacReid(mass, _mass1, _mass2);
  complex<double> mom = PawianQFT::breakupMomQReid(mass, _mass1, _mass2);
  CorrectForChosenSign(mom, rho);
  return rho;
}

complex<double> PhaseSpaceIsobarReid::breakUpMom(const complex<double> mass){
  //  complex<double> q = PawianQFT::breakupMomQReid(mass,_mass1, _mass2);
  complex<double> momReid = PawianQFT::breakupMomQReid(mass, _mass1, _mass2);
  //  complex<double> q = PawianQFT::breakupMomQDefault(mass,_mass1, _mass2);
  complex<double> q = PawianQFT::breakupMomQDefaultAS(mass,_mass1, _mass2);
  CorrectForChosenSign(momReid, q);
  return q;
}

void PhaseSpaceIsobarReid::CorrectForChosenSign(complex<double>& breakUpMom, complex<double>& toChange){
  //assumption relevant imaginary part of the relevant mass is negative
  if(_bumImPartSign > 0 && breakUpMom.imag() > 0.) toChange *= -1;
  if(_bumImPartSign < 0 && !(breakUpMom.imag() > 0.)) toChange *= -1;
  
  //assumption relevant imaginary part of the relevant mass is positive
  // if(_bumImPartSign > 0 && breakUpMom.imag() < 0.) toChange *= -1;
  // if(_bumImPartSign < 0 && !(breakUpMom.imag() < 0.)) toChange *= -1; 
}


