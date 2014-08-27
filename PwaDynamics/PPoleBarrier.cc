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

#include "PwaDynamics/PPoleBarrier.hh"
#include "PwaDynamics/BarrierFactor.hh"

PPoleBarrier::PPoleBarrier(complex<double>& beta, vector<double>& g_i, double mass_0, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs, int orbMom):
  PPole(beta, g_i, mass_0)
  , _phpVecs(phpVecs)
  , _orbMom(orbMom)
{
  _breakUpM0.resize(_phpVecs.size());
  _barrierFactor.resize(_phpVecs.size());
  for(unsigned int i=0; i<_phpVecs.size(); ++i){
    _breakUpM0.at(i)=_phpVecs.at(i)->breakUpMom(_poleMass);
  }
}

PPoleBarrier::~PPoleBarrier(){
}

void PPoleBarrier::evalMatrix(const double mass){

  for (int i=0; i< int(_phpVecs.size()); ++i){
    _barrierFactor.at(i) = BarrierFactor::BlattWeisskopfRatio(_orbMom, _phpVecs.at(i)->breakUpMom(mass), 
							      _breakUpM0.at(i), BarrierFactor::qRDefault);
  }

  double denom=_poleMass*_poleMass-mass*mass;
  if(fabs(denom)<1e-10){
    if(denom<0.) denom=-1e-10; 
    else denom=1e-10;
  }
  for (int i=0; i< int(_g_i.size()); ++i){
    this->operator()(i,0)= (_beta*_g_i.at(i)*_barrierFactor.at(i))/denom;
  }
}

void PPoleBarrier::updatePoleMass (double newPoleMass){
  _poleMass=newPoleMass;
  for(unsigned int i=0; i<_phpVecs.size(); ++i){
    _breakUpM0.at(i)=_phpVecs.at(i)->breakUpMom(_poleMass);
  }
}

