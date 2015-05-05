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

#include "PwaDynamics/KMatrixSlowAdlerCorRel.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"

KMatrixSlowAdlerCorRel::KMatrixSlowAdlerCorRel(vector<std::shared_ptr<KPole> > Kpoles, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs, std::shared_ptr<array_type_2d> fscatPtr, double s0_scat, double sAdler0, double sAdler) :
  KMatrixBase(Kpoles, phpVecs)
  ,_fScatPtr(fscatPtr)
  ,_s0Scat(s0_scat)
  ,_sAdler0(sAdler0)
  ,_sAdler(sAdler)
 {
 }

KMatrixSlowAdlerCorRel::KMatrixSlowAdlerCorRel(int dim) :
  KMatrixBase(dim, dim)
  ,_s0Scat(0.)
  ,_sAdler0(0.)
  ,_sAdler(0.)
{
}

KMatrixSlowAdlerCorRel::~KMatrixSlowAdlerCorRel(){
}


void KMatrixSlowAdlerCorRel::evalMatrix(const double mass){
  evalMatrixTemplate(mass);
}

void KMatrixSlowAdlerCorRel::evalMatrix(const complex<double> mass){
   evalMatrixTemplate(mass);
}

template<typename MassType>
void KMatrixSlowAdlerCorRel::evalMatrixTemplate(const MassType mass){
  for (int i=0; i<NumRows(); ++i){
    for (int j=0; j<NumCols(); ++j){
      this->operator()(i,j)=0.;
    }
  }
  vector<std::shared_ptr<KPole> >::iterator it;
  for (it =_KPoles.begin(); it != _KPoles.end(); ++it){
    (*it)->evalMatrix(mass);
    (*this)+= *(*it);
  }

  MassType sAdler0Denom=mass*mass-_sAdler0;
  double absAdler0Denom=abs(sAdler0Denom);
  if(absAdler0Denom < 1e-10){
     sAdler0Denom *= 1e-10 / absAdler0Denom;
  }

 
  MassType adlerFactor=(1.-_sAdler0)/sAdler0Denom*(mass*mass-_sAdler*0.1349766*0.1349766/2.);

  MassType s0ScatDenom=mass*mass-_s0Scat;
  double absScatDenom=abs(s0ScatDenom);
  if(absScatDenom < 1e-10){
     s0ScatDenom *= 1e-10 / absScatDenom;
  }
  MassType s0ScatFactor=(1.0-_s0Scat)/s0ScatDenom;

  for (int i=0; i<NumRows(); ++i){
    for (int j=0; j<NumCols(); ++j){
      this->operator()(i,j)+=complex<double>((*_fScatPtr)[i][j]*s0ScatFactor);
      this->operator()(i,j)*=adlerFactor;
    }
  }
}

template void KMatrixSlowAdlerCorRel::evalMatrixTemplate(const double mass);
template void KMatrixSlowAdlerCorRel::evalMatrixTemplate(const complex<double> mass);
