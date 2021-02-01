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

#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"

PhaseSpaceIsobar::PhaseSpaceIsobar(double mass1, double mass2):
  AbsPhaseSpace()
  , _mass1(mass1)
  ,_mass2(mass2)
{
  _name="Default";
}

PhaseSpaceIsobar::~PhaseSpaceIsobar(){

}

complex<double> PhaseSpaceIsobar::factor(const double mass, int orbMom){
   return PawianQFT::phaseSpaceFacDefault(mass,_mass1, _mass2);
}

complex<double> PhaseSpaceIsobar::breakUpMom(const double mass){
  return PawianQFT::breakupMomQDefault(mass,_mass1, _mass2);
}

complex<double> PhaseSpaceIsobar::factor(const complex<double> mass, int orbMom){
 //    Calc from the breakup momentum to account for chosen sign
 complex<double> q = PawianQFT::breakupMomQDefault(mass,_mass1, _mass2);
 CorrectForChosenSign(q,q);
 complex<double> rho=2.*q/mass;
 return rho;
}

complex<double> PhaseSpaceIsobar::breakUpMom(const complex<double> mass){
  complex<double> q = PawianQFT::breakupMomQDefault(mass,_mass1, _mass2);
  CorrectForChosenSign(q, q);
  return q;
}

double PhaseSpaceIsobar::thresholdMass(){
  return _mass1+_mass2;
}


