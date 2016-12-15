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

#include "PwaDynamics/Flatte.hh"

#include "Particle/Particle.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"

Flatte::Flatte(std::pair<Particle*, Particle*>& decPair1, std::pair<Particle*, Particle*>& decPair2) :
  _mass11(decPair1.first->mass())
  ,_mass12(decPair1.second->mass())
  ,_mass21(decPair2.first->mass())
  ,_mass22(decPair2.second->mass())
 {
 }

Flatte::Flatte(std::pair <const double, const double>& massPair1, std::pair <const double, const double>& massPair2) :
  _mass11(massPair1.first)
  ,_mass12(massPair1.second)
  ,_mass21(massPair2.first)
  ,_mass22(massPair2.second)
 {
 }
Flatte::~Flatte(){
}

complex<double> Flatte::calcFirstChannel(double currentMass, double mass0, double g1, double g2){
  
  complex<double> i(0.,1.);
  
  //calculate gammas with phase-space factors 
  complex<double> gamma11=g1*PawianQFT::breakupMomQDefault(currentMass, _mass11, _mass12);
  complex<double> gamma22=g2*PawianQFT::breakupMomQDefault(currentMass, _mass21, _mass22);

  complex<double> gammaLow(0.,0.);
  if( (_mass11+_mass12) < (_mass21+_mass22) ) gammaLow=gamma11;
  else gammaLow=gamma22;

  complex<double>  result=mass0*sqrt(gammaLow*gamma11)/( mass0*mass0 - currentMass*currentMass - i * mass0 * (gamma11+gamma22) );
  return result;
}

complex<double> Flatte::calcSecondChannel(double currentMass, double mass0, double g1, double g2){
  
  complex<double> i(0.,1.);
  
  //calculate gammas with phase-space factors 
  complex<double> gamma11=g1*PawianQFT::breakupMomQDefault(currentMass, _mass11, _mass12);
  complex<double> gamma22=g2*PawianQFT::breakupMomQDefault(currentMass, _mass21, _mass22);

  complex<double> gammaLow(0.,0.);
  if( (_mass11+_mass12) < (_mass21+_mass22) ) gammaLow=gamma11;
  else gammaLow=gamma22;

  complex<double>  result=mass0*sqrt(gammaLow*gamma22)/( mass0*mass0 - currentMass*currentMass - i * mass0 * (gamma11+gamma22) );
  return result;
}
