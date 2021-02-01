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
#include "PwaDynamics/PhaseSpaceIsobarReidAngularMomentum.hh"
#include "Utils/PawianConstants.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"

PhaseSpaceIsobarReidAngularMomentum::PhaseSpaceIsobarReidAngularMomentum(double mass1, double mass2):
  PhaseSpaceIsobar(mass1, mass2)
{
  _name="ReidAngularMomentum";
}

PhaseSpaceIsobarReidAngularMomentum::~PhaseSpaceIsobarReidAngularMomentum(){

}

complex<double> PhaseSpaceIsobarReidAngularMomentum::factor(const double mass, int orbMom){
  double expL=2.*orbMom+1.;
  complex<double> cmImag(PawianQFT::ChewMandelstamReid_AngularMomentum(mass*mass, _mass1, _mass2, orbMom).imag(), 0.);
  complex<double> rho=pow(cmImag, expL); 
  return rho;
}

complex<double> PhaseSpaceIsobarReidAngularMomentum::breakUpMom(const double mass){
  return PawianQFT::breakupMomQReid(mass,_mass1, _mass2);
}

complex<double> PhaseSpaceIsobarReidAngularMomentum::factor(const complex<double> mass, int orbMom){
  double expL=2.*orbMom+1.;
  complex<double> cmImag(PawianQFT::ChewMandelstamReid_AngularMomentum(mass*mass, _mass1, _mass2, orbMom).imag(), 0.);
  complex<double> rho=pow(cmImag, expL); 
  return rho;
}

complex<double> PhaseSpaceIsobarReidAngularMomentum::breakUpMom(const complex<double> mass){
  complex<double> momReid = PawianQFT::breakupMomQReid(mass, _mass1, _mass2);
  CorrectForChosenSign(momReid, momReid);
  return momReid;
}

complex<double> PhaseSpaceIsobarReidAngularMomentum::ChewM(const double mass, int orbMom){
  complex<double> massSqrCompl(mass*mass, 1.e-14); // for real s: expansion to s=0 from 1st quadrant
  return PawianQFT::ChewMandelstamReid_AngularMomentum(massSqrCompl, _mass1, _mass2, orbMom);
}

complex<double> PhaseSpaceIsobarReidAngularMomentum::ChewM(const complex<double> mass, int orbMom){
  complex<double> s=mass*mass;
  complex<double> result = PawianQFT::ChewMandelstamReid_AngularMomentum(s, _mass1, _mass2, orbMom);
  complex<double> momReidL = result.imag()*mass/2.;
  CorrectCMForChosenSign(momReidL, result);
  return result;
}

void PhaseSpaceIsobarReidAngularMomentum::CorrectCMForChosenSign(complex<double>& breakUpMom, complex<double>& toChange){
  if((_bumImPartSign > 0 && breakUpMom.imag() < 0) ||
     (_bumImPartSign < 0 && breakUpMom.imag() > 0)){
     toChange = conj(toChange);
  }
}

