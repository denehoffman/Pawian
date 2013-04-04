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

#include "PwaDynamics/FVectorPiPiS.hh"
#include "PwaDynamics/KMatrixPiPiS.hh"
#include "PwaDynamics/PVectorRel.hh"
#include "PwaDynamics/PPole.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include <boost/multi_array.hpp>

FVectorPiPiS::FVectorPiPiS() :
  FVector(5)
 {
   boost::shared_ptr<KMatrixPiPiS> theKMatrix(new KMatrixPiPiS());
   vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs=theKMatrix->phaseSpaceVec();
   vector<boost::shared_ptr<KPole> > kPoles=theKMatrix->kpoles();

   vector<boost::shared_ptr<PPole> > thePpoles;
   complex<double> defaultBeta(1.,0.); 
     
   vector<boost::shared_ptr<KPole> >::iterator it;
   for (it=kPoles.begin(); it!=kPoles.end(); ++it){ 
     std::vector<double> currentGFactors=(*it)->gFactors();
     boost::shared_ptr<PPole> currentPPole(new PPole(defaultBeta, currentGFactors, (*it)->poleMass()));
     thePpoles.push_back(currentPPole);     
   } 

   std::vector<complex <double> > fProdVec;

   for (int i=0; i<int(phpVecs.size()); ++i){
     complex<double> currentVal(1.0,0.);
     fProdVec.push_back(currentVal);
  }

  double s0Prod=-0.0737;
  boost::shared_ptr<PVectorSlowCorRel> thePVector(new PVectorSlowCorRel(thePpoles, phpVecs, fProdVec, s0Prod));   
  _Kmatrix=theKMatrix;
  _Pvector=thePVector;
  _pVectorCor=thePVector;
     
 }


FVectorPiPiS::~FVectorPiPiS(){
}

