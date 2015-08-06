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

#include "PwaDynamics/PVectorKPiSFocus.hh"
#include "PwaDynamics/PPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"
#include "ErrLogger/ErrLogger.hh"

PVectorKPiSFocus::PVectorKPiSFocus(std::shared_ptr<KMatrixKPiSFocus> kMatrix) :
  PVectorRel(kMatrix->phaseSpaceVec())
  ,_kMatrix(kMatrix)
 {
   _aProdVec.resize(2);
   _bProdVec.resize(2);
   _cProdVec.resize(2);
   _phaseProdVec.resize(2);

   complex<double> defaultBeta(1.,0.); 
  
   vector<std::shared_ptr<KPole> > theKPoles=kMatrix->kpoles();
      
   vector<std::shared_ptr<KPole> >::iterator it;
   for (it=theKPoles.begin(); it!=theKPoles.end(); ++it){ 
     std::vector<double> currentGFactors=(*it)->gFactors();
     std::shared_ptr<PPole> currentPPole(new PPole(defaultBeta, currentGFactors, (*it)->poleMass()));
     _Ppoles.push_back(currentPPole);     
   }
 }

PVectorKPiSFocus::~PVectorKPiSFocus(){
}

void PVectorKPiSFocus::evalMatrix(const double mass, Spin OrbMom){

  double s_hat=mass*mass/(_kMatrix->sNorm())-1.;
  
  Matrix< complex<double> > thePVector(NumRows(), 1);
  vector<std::shared_ptr<PPole> >::iterator it;
  for (it =_Ppoles.begin(); it != _Ppoles.end(); ++it){
    (*it)->evalMatrix(mass);
    thePVector += *(*it);
  }

  for (int i=0; i<thePVector.NumRows(); ++i){
    complex<double> expiphiX(cos(_phaseProdVec[i]), sin(_phaseProdVec[i]));  
    this->operator()(i,0)=thePVector(i,0)+(_aProdVec[i]+_bProdVec[i]*s_hat+_cProdVec[i]*s_hat*s_hat)*expiphiX;
  } 

}

void  PVectorKPiSFocus::updateAprod (int i, double aProd){
  if ( i >= int(_aProdVec.size()) ){
    Alert << "aProd parameter with index " << i << " not available!!!" << endmsg;
    exit(0);
  }
  _aProdVec[i]=aProd;
}

void  PVectorKPiSFocus::updateBprod (int i, double bProd){
  if ( i >= int(_bProdVec.size()) ){
    Alert << "bProd parameter with index " << i << " not available!!!" << endmsg;
    exit(0);
  }
  _bProdVec[i]=bProd;
}

void  PVectorKPiSFocus::updateCprod (int i, double cProd){
  if ( i >= int(_cProdVec.size()) ){
    Alert << "cProd parameter with index " << i << " not available!!!" << endmsg;
    exit(0);
  }
  _cProdVec[i]=cProd;
}

void PVectorKPiSFocus::updatePhaseprod (int i, double phaseProd){
  if ( i >= int(_phaseProdVec.size()) ){
    Alert << "phaseProd parameter with index " << i << " not available!!!" << endmsg;
    exit(0);
  }
  _phaseProdVec[i]=phaseProd;
}

