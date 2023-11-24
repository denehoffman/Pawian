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
#include "ErrLogger/ErrLogger.hh"

KPoleBarrier::KPoleBarrier(vector<double>& g_i, double mass_0, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs, int orbMom, bool truncatedBarrier):
  KPole(g_i, mass_0)
  , _phpVecs(phpVecs)
{
  _orbMom=orbMom;
  _woBarrier.resize(_phpVecs.size());
  _truncatedBarrier=truncatedBarrier;
  _breakUpM0.resize(_phpVecs.size());
  _barrierFactor.resize(_phpVecs.size());
  for(unsigned int i=0; i<_phpVecs.size(); ++i){
    // _breakUpM0.at(i)=_phpVecs.at(i)->breakUpMomDefaultAS(mass_0);
    _phpVecs.at(i)->breakUpMom(mass_0);
    _woBarrier.at(i)=false;
    //    if(_phpVecs.at(i)->name()=="ReidAngularMomentum" || _phpVecs.at(i)->name()=="DudekAngularMomentum") _woBarrier.at(i)=true;
  }
}

KPoleBarrier::~KPoleBarrier(){
}



void KPoleBarrier::evalMatrix(const double mass, Spin OrbMom){
  //  evalMatrixTemplate(mass, OrbMom);
 
  for (int i=0; i< int(_phpVecs.size()); ++i){
    if(_woBarrier.at(i)){
      _barrierFactor.at(i) = 1.;
    }
    else if(_truncatedBarrier){
          _barrierFactor.at(i) = BarrierFactor::BlattWeisskopfTensorRatio(OrbMom, _phpVecs.at(i)->breakUpMom(mass), 
                                                                                _breakUpM0.at(i), BarrierFactor::qRDefault);
    }
    else _barrierFactor.at(i) = BarrierFactor::BlattWeisskopfRatio(OrbMom, _phpVecs.at(i)->breakUpMom(mass), 
  							      _breakUpM0.at(i), BarrierFactor::qRDefault);
  }


  double denom=_poleMass*_poleMass-mass*mass;
  if(std::abs(denom) < 1e-10){
    if(denom<0.) denom = -1.E-10;
    else denom = 1.E-10;
  }

  for (int i=0; i< int(_g_i.size()); ++i){
    for (int j=0; j< int(_g_i.size()); ++j){
      this->operator()(i,j)= ( _g_i.at(i)*_g_i.at(j))/denom;
     }
   }
}

//barrier factors not needed for the scan in the complex energy plane
//void KPoleBarrier::evalMatrix(const complex<double> mass, Spin OrbMom){

void KPoleBarrier::updatePoleMass (double newPoleMass){
  _poleMass=newPoleMass;
  _poleMassCompl=complex<double>(_poleMass, 0.);
  for(unsigned int i=0; i<_phpVecs.size(); ++i){
    _breakUpM0.at(i)=_phpVecs.at(i)->breakUpMomDefaultAS(_poleMass);
  }
}

