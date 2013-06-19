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

#include "PwaDynamics/KPoleBarrier.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

KPoleBarrier::KPoleBarrier(vector<double>& g_i, double mass_0, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs, int orbMom):
  KPole(g_i, mass_0)
  , _phpVecs(phpVecs)
  , _orbMom(orbMom)
{
}

KPoleBarrier::~KPoleBarrier(){
}

void KPoleBarrier::evalMatrix(const double mass){

  vector< complex<double> > barrierFactor;

  for (int i=0; i< int(_phpVecs.size()); ++i){
    complex<double> breakUpM=_phpVecs[i]->breakUpMom(mass);
    complex<double> breakUpM0=_phpVecs[i]->breakUpMom(_poleMass);
    barrierFactor.push_back(breakUpM/breakUpM0);
  }

  for (int i=0; i< int(_g_i.size()); ++i){
    for (int j=0; j< int(_g_i.size()); ++j){
      this->operator()(i,j)= ( _g_i[i]*barrierFactor[i]*_g_i[j]*barrierFactor[j])/(_poleMass*_poleMass-mass*mass);
     }
   }
}

