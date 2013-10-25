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

  double adlerFactor=(1.-_sAdler0)/(mass*mass-_sAdler0)*(mass*mass-_sAdler*0.1349766*0.1349766/2.);
  double s0ScatFactor=(1.0-_s0Scat)/(mass*mass-_s0Scat);

  for (int i=0; i<NumRows(); ++i){
    for (int j=0; j<NumCols(); ++j){
      this->operator()(i,j)+=complex<double> ((*_fScatPtr)[i][j]*s0ScatFactor, 0.);
      this->operator()(i,j)*=adlerFactor;
    }
  }
}
