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

#include "PwaDynamics/TMatrixBase.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"


TMatrixBase::TMatrixBase(std::shared_ptr<KMatrixBase> Kmatrix) :
  Matrix< complex<double> >::Matrix(int(Kmatrix->phaseSpaceVec().size()), int(Kmatrix->phaseSpaceVec().size()))
  , _Kmatrix(Kmatrix)
 {
 }


TMatrixBase::~TMatrixBase(){
}

void TMatrixBase::evalMatrix(const double mass){
  _Kmatrix->evalMatrix(mass);
}

void TMatrixBase::evalMatrix(const complex<double> mass){
  _Kmatrix->evalMatrix(mass);
}

void TMatrixBase::SetBumImPartSigns(std::vector<double> signs){
  _Kmatrix->SetBumImPartSigns(signs);
}
