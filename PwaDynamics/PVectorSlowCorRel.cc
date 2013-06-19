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

#include "PwaDynamics/PVectorSlowCorRel.hh"
#include "PwaDynamics/PPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"
#include "ErrLogger/ErrLogger.hh"

PVectorSlowCorRel::PVectorSlowCorRel(vector<std::shared_ptr<PPole> > Ppoles, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs, std::vector< complex<double> >& fProdVec, double s0prod) :
  PVectorRel(Ppoles, phpVecs)
  ,_fProdVec(fProdVec)
  ,_s0prod(s0prod)
 {
 } 

PVectorSlowCorRel::~PVectorSlowCorRel(){
}

void PVectorSlowCorRel::evalMatrix(const double mass){

  Matrix< complex<double> > thePVector(NumRows(), 1);
  vector<std::shared_ptr<PPole> >::iterator it;
  for (it =_Ppoles.begin(); it != _Ppoles.end(); ++it){
    (*it)->evalMatrix(mass);
    thePVector += *(*it);
  }

  for (int i=0; i<thePVector.NumRows(); ++i){  
    this->operator()(i,0)=thePVector(i,0)+ _fProdVec[i]*(1.0-_s0prod)/(mass*mass-_s0prod);
  } 

}

void PVectorSlowCorRel::updateFprod(int i, complex<double> fProd) {
  if ( i >= int(_fProdVec.size()) ){
    Alert << "no fProd parameter with index " << i << " available!!!" << endmsg;
    exit(0);
  }
  _fProdVec[i]=fProd;
}
