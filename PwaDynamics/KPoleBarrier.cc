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
#include "PwaDynamics/BarrierFactor.hh"

KPoleBarrier::KPoleBarrier(vector<double>& g_i, double mass_0, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs, int orbMom, bool truncatedBarrier):
  KPole(g_i, mass_0)
  , _phpVecs(phpVecs)
  , _orbMom(orbMom)
  ,_truncatedBarrier(truncatedBarrier)
{
  _breakUpM0.resize(_phpVecs.size());
  _barrierFactor.resize(_phpVecs.size());
  for(unsigned int i=0; i<_phpVecs.size(); ++i){
    _breakUpM0.at(i)=_phpVecs.at(i)->breakUpMom(mass_0);
  }
}

KPoleBarrier::~KPoleBarrier(){
}



void KPoleBarrier::evalMatrix(const double mass, Spin OrbMom){
  evalMatrixTemplate(mass, OrbMom);
}

void KPoleBarrier::evalMatrix(const complex<double> mass, Spin OrbMom){
  evalMatrixTemplate(mass, OrbMom);
}
 
template<typename MassType>
void KPoleBarrier::evalMatrixTemplate(const MassType mass, Spin OrbMom){

  for (int i=0; i< int(_phpVecs.size()); ++i){
    if(_truncatedBarrier){
      _barrierFactor.at(i) = BarrierFactor::BlattWeisskopfTensorRatio(OrbMom, _phpVecs.at(i)->breakUpMom(mass), 
								      _breakUpM0.at(i), BarrierFactor::qRDefault);
    }
    else{
      _barrierFactor.at(i) = BarrierFactor::BlattWeisskopfRatio(OrbMom, _phpVecs.at(i)->breakUpMom(mass), 
							       _breakUpM0.at(i), BarrierFactor::qRDefault);
    }
  }

  MassType denom=_poleMass*_poleMass-mass*mass;
  double absDenom=abs(denom);
  if(absDenom < 1e-10){
     denom *= 1E-10/absDenom;
  }

  for (int i=0; i< int(_g_i.size()); ++i){
    for (int j=0; j< int(_g_i.size()); ++j){
      this->operator()(i,j)= ( _g_i.at(i)*_barrierFactor.at(i)*_g_i.at(j)*_barrierFactor.at(j))/denom;
     }
   }
}

template void KPoleBarrier::evalMatrixTemplate(const double mass, Spin OrbMom);
template void KPoleBarrier::evalMatrixTemplate(const complex<double> mass, Spin OrbMom);

void KPoleBarrier::updatePoleMass (double newPoleMass){
  _poleMass=newPoleMass;
  for(unsigned int i=0; i<_phpVecs.size(); ++i){
    _breakUpM0.at(i)=_phpVecs.at(i)->breakUpMom(_poleMass);
  }
}

