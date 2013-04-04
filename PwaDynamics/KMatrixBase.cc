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

#include "PwaDynamics/KMatrixBase.hh"
#include "PwaDynamics/KPole.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"

KMatrixBase::KMatrixBase(vector<boost::shared_ptr<KPole> > Kpoles, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs) :
  Matrix< complex<double> >::Matrix(int(phpVecs.size()), int(phpVecs.size()))
  ,_KPoles(Kpoles)
  ,_phpVecs(phpVecs)
 {
 }

KMatrixBase::KMatrixBase(vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs, int numCols, int numRows) :
  Matrix< complex<double> >::Matrix(numCols, numRows)
  ,_phpVecs(phpVecs)
 {
 }
KMatrixBase::KMatrixBase(int numCols, int numRows) :
  Matrix< complex<double> >::Matrix(numCols, numRows)
 {
 }

KMatrixBase::~KMatrixBase(){
}

