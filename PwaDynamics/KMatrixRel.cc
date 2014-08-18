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

#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"

KMatrixRel::KMatrixRel(vector<std::shared_ptr<KPole> > Kpoles, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs) :
  KMatrixBase(Kpoles, phpVecs)
 {
 }


KMatrixRel::~KMatrixRel(){
}

void KMatrixRel::evalMatrix(const double mass){
   evalMatrixTemplate(mass);
}

void KMatrixRel::evalMatrix(const complex<double> mass){
   evalMatrixTemplate(mass);
}

template<typename MassType>
void KMatrixRel::evalMatrixTemplate(const MassType mass){

  Matrix< complex<double> > theKMatrix(NumRows(), NumRows());
  vector<std::shared_ptr<KPole> >::iterator it;
  for (it =_KPoles.begin(); it != _KPoles.end(); ++it){
    (*it)->evalMatrix(mass);
    theKMatrix += *(*it);
  }

  for (int i=0; i<theKMatrix.NumRows(); ++i){
    for (int j=0; j<theKMatrix.NumCols(); ++j){
      this->operator()(i,j)=theKMatrix(i,j);
    }
  }
}

template void KMatrixRel::evalMatrixTemplate(const double mass);
template void KMatrixRel::evalMatrixTemplate(const complex<double> mass);
