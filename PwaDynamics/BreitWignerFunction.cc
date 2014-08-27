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

#include "PwaDynamics/BreitWignerFunction.hh"
#include "PwaDynamics/BarrierFactor.hh"

complex<double>  BreitWignerFunction::NonRel(double currentMass,double mass0, double width){
  complex<double> i(0.,1.);
  return mass0*width/(mass0*mass0 - currentMass*currentMass- i*mass0*width);
}

complex<double> BreitWignerFunction::Rel(double currentMass,double mass0, double width, double massA, double massB){
  complex<double> i(0.,1.);
  complex<double> rho0_m=phaseSpaceFac(mass0, massA, massB);
  return mass0*width/(mass0*mass0 - currentMass*currentMass- i*rho0_m*mass0*width);
}

complex<double>  BreitWignerFunction::BlattWRel(int orbMom, double currentMass,double mass0, double width, double massA, double massB){
  complex<double> i(0.,1.);
  complex<double> rho0=phaseSpaceFac(mass0, massA, massB);
  complex<double> rho=phaseSpaceFac(currentMass, massA, massB);
  complex<double> momQ0=breakupMomQ(mass0, massA, massB);
  complex<double> momQ=breakupMomQ(currentMass, massA, massB);

  std::complex<double> bwkFactor=BarrierFactor::BlattWeisskopfRatio(orbMom, momQ, momQ0, BarrierFactor::qRDefault);


  return mass0*width*bwkFactor/(mass0*mass0 - currentMass*currentMass- i*(rho/rho0)*mass0*width*bwkFactor*bwkFactor);
  //  return mass0*width*bwkFactor/(mass0*mass0 - currentMass*currentMass- i*mass0*width*bwkFactor*bwkFactor);
  //  return mass0*width/(mass0*mass0 - currentMass*currentMass- i*mass0*width*(rho/rho0)); 
}
