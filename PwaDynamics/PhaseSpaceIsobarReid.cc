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
#include "Utils/PawianConstants.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"

PhaseSpaceIsobarReid::PhaseSpaceIsobarReid(double mass1, double mass2):
  PhaseSpaceIsobar(mass1, mass2)
{
}

PhaseSpaceIsobarReid::~PhaseSpaceIsobarReid(){

}

complex<double> PhaseSpaceIsobarReid::factor(const double mass){
   return ChewM(mass).imag();
}

complex<double> PhaseSpaceIsobarReid::breakUpMom(const double mass){
  return PawianQFT::breakupMomQReid(mass,_mass1, _mass2);
}

complex<double> PhaseSpaceIsobarReid::factor(const complex<double> mass){
  complex<double> rho = ChewM(mass).imag(); 
  return rho;
}

complex<double> PhaseSpaceIsobarReid::breakUpMom(const complex<double> mass){
  complex<double> momReid = PawianQFT::breakupMomQReid(mass, _mass1, _mass2);
  CorrectForChosenSign(momReid, momReid);
  return momReid;
}

complex<double> PhaseSpaceIsobarReid::ChewM(const double mass){
  complex<double> massSqrCompl(mass*mass, 1.e-10); // for real s: expansion to s=0 from 1st quadrant
  return PawianQFT::ChewMandelstamReid(massSqrCompl, _mass1, _mass2);  
}

complex<double> PhaseSpaceIsobarReid::ChewM(const complex<double> mass){
  complex<double> s=mass*mass;
  complex<double> result = PawianQFT::ChewMandelstamReid(s, _mass1, _mass2);
  complex<double> momReid = PawianQFT::breakupMomQReid(mass, _mass1, _mass2);
  CorrectForChosenSign(momReid, result);
  return result;
}

