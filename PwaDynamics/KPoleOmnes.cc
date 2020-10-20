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

#include "PwaDynamics/KPoleOmnes.hh"
#include "PwaDynamics/KPole.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

KPoleOmnes::KPoleOmnes(vector<double>& g_i, double mass_0):
  KPole(g_i, mass_0)
{
}

KPoleOmnes::KPoleOmnes(vector<double>& g_i, double mass_0, int numRow, int numCol):
  KPole(g_i, mass_0, numRow, numCol)
{
}

KPoleOmnes::~KPoleOmnes(){
}

void KPoleOmnes::evalMatrix(const double mass, Spin OrbMom){
  //   evalMatrixTemplate(mass);
  double denom=mass*mass*(mass*mass-_poleMass*_poleMass);

  if( std::abs(denom) < 1.e-10){
      if(denom<0.) denom = -1.e-10;
      else denom = 1.e-10;
  }

  double s_denom=_poleMass*_poleMass/denom;

  for (int i=0; i< int(_g_i.size()); ++i){
    for (int j=0; j< int(_g_i.size()); ++j){
      this->operator()(i,j)= ( _g_i[i]*s_denom*_g_i[j]);
     }
  }
}

void KPoleOmnes::evalMatrix(const complex<double> mass, Spin OrbMom){
  //   evalMatrixTemplate(mass);
  complex<double> denom=mass*mass*(mass*mass-_poleMassCompl*_poleMassCompl);

  if( std::abs(denom) < 1.e-10){
     denom = complex<double>(1.E-10, 0.);
  }

  complex<double> s_denom=_poleMassCompl*_poleMassCompl/denom;
  for (int i=0; i< int(_g_i.size()); ++i){
    for (int j=0; j< int(_g_i.size()); ++j){
      this->operator()(i,j)= ( _g_i[i]*s_denom*_g_i[j]);
     }
  }
}


