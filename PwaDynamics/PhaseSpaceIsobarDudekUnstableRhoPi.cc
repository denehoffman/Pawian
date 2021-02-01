//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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
#include "PwaDynamics/PhaseSpaceIsobarDudekUnstableRhoPi.hh"
#include "Utils/PawianConstants.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"

PhaseSpaceIsobarDudekUnstableRhoPi::PhaseSpaceIsobarDudekUnstableRhoPi(double mass1, double mass2):
  PhaseSpaceIsobar(mass1, mass2)
{
  _name="DudekUnstableRhoPi";
}

PhaseSpaceIsobarDudekUnstableRhoPi::~PhaseSpaceIsobarDudekUnstableRhoPi(){

}

complex<double> PhaseSpaceIsobarDudekUnstableRhoPi::factor(const double mass, int orbMom){
  complex<double> rho(-ChewM(mass).imag(), 0.); 
  return rho;
}

complex<double> PhaseSpaceIsobarDudekUnstableRhoPi::breakUpMom(const double mass){
  return PawianQFT::breakupMomQDudek(mass,_mass1, _mass2);
}

complex<double> PhaseSpaceIsobarDudekUnstableRhoPi::factor(const complex<double> mass, int orbMom){
  complex<double> rho(-ChewM(mass).imag(), 0.); 
  return rho;
}

complex<double> PhaseSpaceIsobarDudekUnstableRhoPi::breakUpMom(const complex<double> mass){
  complex<double> momDudek = PawianQFT::breakupMomQDudek(mass, _mass1, _mass2);
  CorrectForChosenSign(momDudek, momDudek);
  return momDudek;
}

complex<double> PhaseSpaceIsobarDudekUnstableRhoPi::ChewM(const double mass, int orbMom){
  double s=mass*mass;
  return PawianQFT::ChewMandelstamDudekUnstableRhoPi(s, _mass1, _mass2);  
}

complex<double> PhaseSpaceIsobarDudekUnstableRhoPi::ChewM(const complex<double> mass, int orbMom){
  complex<double> result(0.,0.); //dummy
  //  complex<double> s=mass*mass;
  //  complex<double> result = PawianQFT::ChewMandelstamDudek(s, _mass1, _mass2);
  //  complex<double> momDudek = -result.imag()*mass/2.;
  //  CorrectCMForChosenSign(momDudek, result);
  return result;
}

void PhaseSpaceIsobarDudekUnstableRhoPi::CorrectCMForChosenSign(complex<double>& breakUpMom, complex<double>& toChange){
  if((_bumImPartSign > 0 && breakUpMom.imag() < 0) ||
     (_bumImPartSign < 0 && breakUpMom.imag() > 0)){
     toChange = conj(toChange);
  }
}

