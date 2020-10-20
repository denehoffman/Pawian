//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#include "PwaDynamics/PPoleOmnes.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

PPoleOmnes::PPoleOmnes(complex<double>& beta, vector<double>& g_i, double mass_0):
  PPole(beta, g_i, mass_0)
{
}

PPoleOmnes::~PPoleOmnes(){
}

void PPoleOmnes::evalMatrix(const double mass, Spin OrbMom){
  double denom=_poleMass*_poleMass-mass*mass;
  if(fabs(denom)<1.e-10){
    if(denom<0.) denom=-1.e-10; 
    else denom=1.e-10;
  }

  double s_denom=mass*mass/denom;
  
  for (int i=0; i< int(_g_i.size()); ++i){
    this->operator()(i,0)= (_beta*s_denom*_g_i[i]);
  }
}

void PPoleOmnes::evalMatrix(const complex<double> mass, Spin OrbMom){
  complex<double> denom=_poleMass*_poleMass-mass*mass;
  complex<double> s_denom= mass*mass/denom;
  
  for (int i=0; i< int(_g_i.size()); ++i){
    this->operator()(i,0)= _beta*s_denom*_g_i[i];
  }
}

