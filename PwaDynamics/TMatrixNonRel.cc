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

#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/KMatrixNonRel.hh"
#include "qft++/matrix/IdentityMatrix.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"


TMatrixNonRel::TMatrixNonRel(std::shared_ptr<KMatrixNonRel> Kmatrix) :
  TMatrixBase(Kmatrix)
 {
 }


TMatrixNonRel::~TMatrixNonRel(){
}

void TMatrixNonRel::evalMatrix(const double mass){
   evalMatrixTemplate(mass);
}

void TMatrixNonRel::evalMatrix(const complex<double> mass){
   evalMatrixTemplate(mass);
}

template<typename MassType>
void TMatrixNonRel::evalMatrixTemplate(const MassType mass){
  TMatrixBase::evalMatrix(mass);

  complex<double> imagCompl(0.,1.);  
  IdentityMatrix< complex<double> > theIdMatrix(NumRows());
  Matrix< complex< double > > tmpDenomMatrCompl = theIdMatrix-imagCompl*(*_Kmatrix);
  
  Matrix< complex< double > > tmpDenomMatrInv=tmpDenomMatrCompl;
  tmpDenomMatrInv.invert();
  
  Matrix< complex <double> > currentTMatr=tmpDenomMatrInv*(*_Kmatrix);
  
  for (int i=0; i<currentTMatr.NumRows(); ++i){
    for (int j=0; j<currentTMatr.NumCols(); ++j){
      this->operator()(i,j)=currentTMatr(i,j);
    }
  }
}

template void TMatrixNonRel::evalMatrixTemplate(const double mass);
template void TMatrixNonRel::evalMatrixTemplate(const complex<double> mass);
