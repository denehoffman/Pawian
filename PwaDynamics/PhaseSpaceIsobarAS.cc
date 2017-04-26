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

#include "PwaDynamics/PhaseSpaceIsobarAS.hh"
#include "qft++Extension/PawianUtils.hh"

PhaseSpaceIsobarAS::PhaseSpaceIsobarAS(double mass1, double mass2):
  AbsPhaseSpace()
  , _mass1(mass1)
  ,_mass2(mass2)
{
}

PhaseSpaceIsobarAS::~PhaseSpaceIsobarAS(){

}

complex<double> PhaseSpaceIsobarAS::factor(const double mass){
  return PawianQFT::phaseSpaceFacAS(mass,_mass1, _mass2);
}

complex<double> PhaseSpaceIsobarAS::breakUpMom(const double mass){
  return PawianQFT::breakupMomQAS(mass,_mass1, _mass2);
}

complex<double> PhaseSpaceIsobarAS::factor(const complex<double> mass){
  // Calc from the breakup momentum to account for chosen sign
  complex<double> q = PawianQFT::breakupMomQAS(mass,_mass1, _mass2);
  CorrectForChosenSign(q,q);
  return q * 2. / mass;
}

complex<double> PhaseSpaceIsobarAS::breakUpMom(const complex<double> mass){
   
  complex<double> q = PawianQFT::breakupMomQAS(mass,_mass1, _mass2);
  CorrectForChosenSign(q, q);
  return q;
}

double PhaseSpaceIsobarAS::thresholdMass(){
  return _mass1+_mass2;
}

