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

#include "PwaDynamics/FVector.hh"
#include "ErrLogger/ErrLogger.hh"
//#include "qft++/relativistic-quantum-mechanics/Utils.hh"


FVector::FVector(std::shared_ptr<KMatrixBase> Kmatrix, std::shared_ptr<PVectorRel> Pvector) :
  Matrix< complex<double> >::Matrix(Pvector->NumRows(), Pvector->NumCols())
  , _Kmatrix(Kmatrix)
  ,_Pvector(Pvector)
  ,_imagCompl(0.,1.)
  ,_idMatrix(NumRows())
  ,_CMMatrix(NumRows(),NumRows())
  ,_phpVec(Kmatrix->phaseSpaceVec())
 {
   for (int i=0; i<NumRows(); ++i){
     for (int j=0; j<NumRows(); ++j){
       _CMMatrix(i,j)=complex<double> (0.,0.);
     }
   }
 }

FVector::FVector(int numRows) :
  Matrix< complex<double> >::Matrix(numRows, 1)
  ,_imagCompl(0.,1.)
  ,_idMatrix(NumRows())
  ,_CMMatrix(NumRows(),NumRows())
 {
   for (int i=0; i<NumRows(); ++i){
     for (int j=0; j<NumRows(); ++j){
       _CMMatrix(i,j)=complex<double> (0.,0.);
     }
   }
 }

FVector::~FVector(){
}

void FVector::evalMatrix(const double mass, Spin OrbMom){

  _Kmatrix->evalMatrix(mass, OrbMom);
  _Pvector->evalMatrix(mass, OrbMom);

 // for (int i=0; i<NumRows(); ++i) _rhoMatrix(i,i) = _phpVec[i]->factor(mass);
  for (int i=0; i<NumRows(); ++i) _CMMatrix(i,i) = _phpVec[i]->ChewM(mass);

  Matrix< complex< double > > denomMatrComplInv = _idMatrix-(*_Kmatrix)*_CMMatrix;
  denomMatrComplInv.invert();
  
  Matrix< complex <double> > currentTMatr=denomMatrComplInv*(*_Pvector);

  for (int i=0; i<currentTMatr.NumRows(); ++i){
    this->operator()(i,0)=currentTMatr(i,0);
  }

}

complex<double> FVector::evalProjMatrix(const double mass, int index, Spin OrbMom){
  _Kmatrix->evalMatrix(mass, OrbMom);
  _Pvector->evalMatrix(mass, OrbMom);

  for (int i=0; i<NumRows(); ++i) _CMMatrix(i,i) = _phpVec[i]->ChewM(mass);
 
  Matrix< complex< double > > denomMatrInv=_idMatrix-(*_Kmatrix)*_CMMatrix;
  denomMatrInv.invert();
  
  complex <double> result(0.,0.);
  for(int i=0;i<NumRows(); ++i){
    result+=denomMatrInv(index,i)*(*_Pvector)(i,0);
  }

  if(result.real() != result.real()){
    Alert << "result:\t" << result << endmsg;
    Alert << "mass:\t" << mass << "\tindex: " << index <<"\tOrbMom: " << OrbMom << endmsg;
    for(int i=0;i<NumRows(); ++i){
      Alert <<"(*_Pvector)(" << i << ",0): " << (*_Pvector)(i,0) << endmsg;
      Alert <<"denomMatrInv(index," << i << "): " << denomMatrInv(index,i) << endmsg;
      Alert << "_Kmatrix(" << i << "," << i << "): " << (*_Kmatrix)(i,i) << endmsg;
      Alert << "_CMMatrix(" << i << "," << i << "): " << _CMMatrix(i,i) << endmsg;
    }
    exit(0);
  }
  return result; 
}

