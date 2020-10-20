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

#include "PwaDynamics/KMatrixOmnes.hh"
#include "PwaDynamics/KPoleOmnes.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"

KMatrixOmnes::KMatrixOmnes(vector<std::shared_ptr<KPole> > KpolesOmnes, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs) : KMatrixNonRel(KpolesOmnes, phpVecs)
 {
 }


KMatrixOmnes::~KMatrixOmnes(){
}

void KMatrixOmnes::evalMatrix(const double mass, Spin OrbMom){
  return evalMatrixTemplate(mass, OrbMom);
}

void KMatrixOmnes::evalMatrix(const complex<double> mass, Spin OrbMom){
   return evalMatrixTemplate(mass, OrbMom);
}

template<typename MassType>
void KMatrixOmnes::evalMatrixTemplate(const MassType mass, Spin OrbMom){

  Matrix< complex<double> > theKMatrix(NumRows(), NumRows());
  vector<std::shared_ptr<KPole> >::iterator it;
  for (it =_KPoles.begin(); it != _KPoles.end(); ++it){
    (*it)->evalMatrix(mass, OrbMom);
    theKMatrix += *(*it);
  }

  for (int i=0; i<theKMatrix.NumRows(); ++i){
    for (int j=0; j<theKMatrix.NumCols(); ++j){
      this->operator()(i,j)=theKMatrix(i,j);
    }
  }
}

template void KMatrixOmnes::evalMatrixTemplate(const double mass, Spin OrbMom);
template void KMatrixOmnes::evalMatrixTemplate(const complex<double> mass, Spin OrbMom);
