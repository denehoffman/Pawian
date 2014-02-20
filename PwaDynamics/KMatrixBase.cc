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
#include "ErrLogger/ErrLogger.hh"

KMatrixBase::KMatrixBase(vector<std::shared_ptr<KPole> > Kpoles, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs) :
  Matrix< complex<double> >::Matrix(int(phpVecs.size()), int(phpVecs.size()))
  ,_KPoles(Kpoles)
  ,_phpVecs(phpVecs)
  ,_orderBg(0)
  ,_s0Adler(0.)
  ,_snormAdler(1.)
 {
 }

KMatrixBase::KMatrixBase(vector<std::shared_ptr<AbsPhaseSpace> > phpVecs, int numCols, int numRows) :
  Matrix< complex<double> >::Matrix(numCols, numRows)
  ,_phpVecs(phpVecs)
  ,_orderBg(0)
  ,_s0Adler(0.)
  ,_snormAdler(1.)
 {
 }

KMatrixBase::KMatrixBase(int numCols, int numRows) :
  Matrix< complex<double> >::Matrix(numCols, numRows)
  ,_orderBg(0)
 {
 }

KMatrixBase::~KMatrixBase(){
}

void KMatrixBase::updateBgTerms(unsigned int order, unsigned int row,  unsigned int column, double theVal){
  if(row>column){
    Alert << "K-matrix is symmetric; for update row<=column" << endmsg;
    exit(0);
  }
  if(order>_orderBg){
    Alert << "background parameter for order " << order << " not available!!!" << endmsg;
    exit(0);
  }
  if(int(row) >= NumRows() || int(column) >= NumCols()){
    Alert << "row " << row << " or column " << column 
	  << " >= NumRows " << NumRows() << " or  >= NumCols " << NumCols() << endmsg;
    exit(0);
  }

  _bgTerms.at(order).at(row).at(column)=theVal;
  _bgTerms.at(order).at(column).at(row)=theVal;
}

