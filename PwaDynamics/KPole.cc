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

#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

KPole::KPole(vector<double>& g_i, double mass_0):
  Matrix< complex<double> >::Matrix(int(g_i.size()), int(g_i.size()))
  , _g_i(g_i)
  , _poleMass(mass_0)
{
}

KPole::KPole(vector<double>& g_i, double mass_0, int numRow, int numCol):
  Matrix< complex<double> >::Matrix(numRow, numCol)
  , _g_i(g_i)
  , _poleMass(mass_0)
{
}

KPole::~KPole(){
}

void KPole::evalMatrix(const double mass, Spin OrbMom){
   evalMatrixTemplate(mass);
}

void KPole::evalMatrix(const complex<double> mass, Spin OrbMom){
   evalMatrixTemplate(mass);
}

template<typename MassType>
void KPole::evalMatrixTemplate(const MassType mass, Spin OrbMom){
  MassType denom=_poleMass*_poleMass-mass*mass;

  double absDenom = abs(denom);
  if(absDenom < 1e-10){
     denom *= 1e-10 / absDenom;
  }

  for (int i=0; i< int(_g_i.size()); ++i){
    for (int j=0; j< int(_g_i.size()); ++j){
      this->operator()(i,j)= ( _g_i[i]*_g_i[j])/denom;
     }
  }
}

template void KPole::evalMatrixTemplate(const double mass, Spin OrbMom);
template void KPole::evalMatrixTemplate(const complex<double> mass, Spin OrbMom);
