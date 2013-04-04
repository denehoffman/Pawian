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

#include "PwaDynamics/PPole.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

PPole::PPole(complex<double>& beta, vector<double>& g_i, double mass_0):
  KPole(g_i, mass_0, g_i.size(), 1)
  ,_beta(beta)
{
}

PPole::~PPole(){
}

void PPole::evalMatrix(const double mass){

  for (int i=0; i< int(_g_i.size()); ++i){
    this->operator()(i,0)= (_beta*_g_i[i])/(_poleMass*_poleMass-mass*mass);
  }
}

