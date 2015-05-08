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

#include "PwaDynamics/RadMultipoleFormFactor.hh"

// Considering a decay of the form: A -> gamma B
complex<double> RadMultipoleFormFactor::PureM1(double massA, double currentMassB, double massB0, double Egamma){
  complex<double> i(0.,1.);
  double E0 = (massA*massA - massB0*massB0)/(2.*massA);
  double fD = (E0*E0)/(Egamma*E0 + ((Egamma-E0)*(Egamma-E0)));
  complex<double> M1(Egamma*Egamma*Egamma*fD,0.0);

  return sqrt(M1);
}
