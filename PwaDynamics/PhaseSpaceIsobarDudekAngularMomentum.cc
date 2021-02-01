//************************************************************************//
//									  //
//  Copyright 2021 Bertram Kopf (bertram@ep1.rub.de)			  //
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
#include "PwaDynamics/PhaseSpaceIsobarDudekAngularMomentum.hh"
#include "Utils/PawianConstants.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"

PhaseSpaceIsobarDudekAngularMomentum::PhaseSpaceIsobarDudekAngularMomentum(double mass1, double mass2):
  PhaseSpaceIsobarDudek(mass1, mass2)
{
  _name="DudekAngularMomentum";
}

PhaseSpaceIsobarDudekAngularMomentum::~PhaseSpaceIsobarDudekAngularMomentum(){

}

complex<double> PhaseSpaceIsobarDudekAngularMomentum::factor(const double mass, int orbMom){
  double expL=2.*orbMom+1.;
  complex<double> cmImag(-PawianQFT::ChewMandelstamReid(mass*mass, _mass1, _mass2).imag(), 0.);
  complex<double> rho=std::pow(cmImag, expL); 
  return rho;
}

complex<double> PhaseSpaceIsobarDudekAngularMomentum::factor(const complex<double> mass, int orbMom){
  double expL=2.*orbMom+1.;
  complex<double> cmImag(-PawianQFT::ChewMandelstamReid(mass*mass, _mass1, _mass2).imag(), 0);
  complex<double> rho=std::pow(cmImag, expL); 
  return rho;
}

complex<double> PhaseSpaceIsobarDudekAngularMomentum::ChewM(const double mass, int orbMom){
  complex<double> massSqrCompl(mass*mass, 1.e-14); // for real s: expansion to s=0 from 1st quadrant
  return PawianQFT::ChewMandelstamDudek_AngularMomentum(massSqrCompl, _mass1, _mass2, orbMom);  
}

complex<double> PhaseSpaceIsobarDudekAngularMomentum::ChewM(const complex<double> mass, int orbMom){
  complex<double> s=mass*mass;
  complex<double> result = PawianQFT::ChewMandelstamDudek_AngularMomentum(s, _mass1, _mass2, orbMom);
  complex<double> momDudek = -result.imag()*mass/2.;
  CorrectCMForChosenSign(momDudek, result);
  return result;
}

