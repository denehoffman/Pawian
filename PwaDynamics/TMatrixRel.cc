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

#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/KMatrixBase.hh"
#include "Utils/PawianConstants.hh"
#include "qft++/matrix/IdentityMatrix.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"


TMatrixRel::TMatrixRel(std::shared_ptr<KMatrixBase> Kmatrix) :
  TMatrixBase(Kmatrix)
 {
 }


TMatrixRel::~TMatrixRel(){
}


void TMatrixRel::evalMatrix(const double mass){
   evalMatrixTemplate(mass);
}

void TMatrixRel::evalMatrix(const complex<double> mass){
   evalMatrixTemplate(mass);
}

template<typename MassType>
void TMatrixRel::evalMatrixTemplate(const MassType mass){

   TMatrixBase::evalMatrix(mass);

   vector<std::shared_ptr<AbsPhaseSpace> > phpVec=_Kmatrix->phaseSpaceVec();
   //   Matrix< complex<double> > theRhoMatrix(NumRows(),NumRows());
   Matrix< complex<double> > ChewMMatrix(NumRows(),NumRows());

   for (int i=0; i<NumRows(); ++i){
      for (int j=0; j<NumRows(); ++j){
	 ChewMMatrix(i,j)=complex<double> (0.,0.);
	 if (i==j){
	   //	    theRhoMatrix(i,j) = phpVec[j]->factor(mass);
	   ChewMMatrix(i,j) = phpVec[j]->ChewM(mass);
	 }
      }
   }

   IdentityMatrix< complex<double> > theIdMatrix(NumRows());
   Matrix< complex< double > > tmpDenomMatrInv = theIdMatrix-(*_Kmatrix)*ChewMMatrix;

   tmpDenomMatrInv.invert();

   Matrix< complex <double> > currentTMatr=tmpDenomMatrInv*(*_Kmatrix);

   for (int i=0; i<currentTMatr.NumRows(); ++i){
      for (int j=0; j<currentTMatr.NumCols(); ++j){
	 this->operator()(i,j)=currentTMatr(i,j);
      }
   }
}



void TMatrixRel::evalLorenzInvMatrix(const complex<double> mass){
   vector<std::shared_ptr<AbsPhaseSpace> > phpVec=_Kmatrix->phaseSpaceVec();
   evalMatrix(mass);
   for (int i=0; i<NumRows(); ++i){
      for (int j=0; j<NumCols(); ++j){
	this->operator()(i,j)=std::sqrt(phpVec.at(i)->factor(mass))*this->operator()(i,j)*std::sqrt(phpVec.at(j)->factor(mass));
      }
   }
}

template void TMatrixRel::evalMatrixTemplate(const double);
template void TMatrixRel::evalMatrixTemplate(const complex<double>);
