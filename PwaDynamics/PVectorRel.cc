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

#include "PwaDynamics/PVectorRel.hh"
#include "PwaDynamics/PPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"
#include "ErrLogger/ErrLogger.hh"

PVectorRel::PVectorRel(vector<std::shared_ptr<PPole> > Ppoles, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs) :
  KMatrixBase(phpVecs, int(phpVecs.size()), 1)
  ,_Ppoles(Ppoles)
  ,_nOfChannels(phpVecs.size())
 {
 }

PVectorRel::PVectorRel(vector<std::shared_ptr<AbsPhaseSpace> > phpVecs) :
  KMatrixBase(phpVecs, int(phpVecs.size()), 1)
  ,_nOfChannels(phpVecs.size())
 {
 }

PVectorRel::~PVectorRel(){
}

void PVectorRel::evalMatrix(const double mass, Spin OrbMom){

  Matrix< complex<double> > thePVector(NumRows(), 1);
  vector<std::shared_ptr<PPole> >::iterator it;
  for (it =_Ppoles.begin(); it != _Ppoles.end(); ++it){
    (*it)->evalMatrix(mass, OrbMom);
    thePVector += *(*it);
  }

  for (int i=0; i<thePVector.NumRows(); ++i){
      this->operator()(i,0)=thePVector(i,0);
  }



}

void PVectorRel::updateBeta(int i, complex<double> beta){
  if ( i<0 || i>=int(_Ppoles.size()) ){
    Alert << "no pole with index " << i << " available!!!" << endmsg;
    exit(0);
  }
  _Ppoles[i]->updateBeta(beta);
}

void PVectorRel::updatePoleMass(int i, double mass){
  if ( i<0 || i>=int(_Ppoles.size()) ){
    Alert << "no pole with index " << i << " available!!!" << endmsg;
    exit(0);
  }
  _Ppoles[i]->updatePoleMass(mass);
}

void PVectorRel::updategFactors(int i, vector<double>& newg_i){
  if ( i<0 || i>=int(_Ppoles.size()) ){
    Alert << "no pole with index " << i << " available!!!" << endmsg;
    exit(0);
  }
  _Ppoles[i]->updategFactors(newg_i);
}

void PVectorRel::printElements(){
  std::cout << "\n***P-vector*****" << std::endl;
  for(int i=0; i< NumRows(); ++i){
      std::cout << i << " = " << this->operator()(i,0) << std::endl;
  }
}
