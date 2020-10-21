//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#include "PwaDynamics/FVectorOmnes.hh"
#include "PwaDynamics/KMatrixOmnes.hh"
#include "PwaDynamics/PVectorBgOmnes.hh"
#include "PwaDynamics/PhaseSpaceList.hh"
#include "ErrLogger/ErrLogger.hh"


//#include "qft++/relativistic-quantum-mechanics/Utils.hh"


FVectorOmnes::FVectorOmnes(std::shared_ptr<KMatrixOmnes> Kmatrix, std::shared_ptr<PVectorBgOmnes> pVecBgOmnes, std::shared_ptr<PhaseSpaceList> phpList) :
  FVector(Kmatrix, pVecBgOmnes)
  , _phpList(phpList)
  , _OmnesComplDMatrix(NumRows(),NumRows())
  , _SelfEnergyComplDMatrix(NumRows(),NumRows())
 {
 }

//FVectorOmnes::FVectorOmnes(int numRows) :
//  FVector(numRows)
// {
// }

FVectorOmnes::~FVectorOmnes(){
}

void FVectorOmnes::evalMatrix(const double mass, Spin OrbMom){

  _Kmatrix->evalMatrix(mass, OrbMom);
  _Pvector->evalMatrix(mass, OrbMom);

  for (int i=0; i<NumRows(); ++i){
    for (int j=0; j<NumRows(); ++j){
      _OmnesComplDMatrix(i,j) = _phpList->omnesMatrix()[i][j]->eval(mass*mass);
      _SelfEnergyComplDMatrix(i,j) = _phpList->selfEnergyMatrix()[i][j]->eval(mass*mass);
    }
  }
  Matrix< complex< double > > denomMatrComplInv = _idMatrix+(*_Kmatrix)*_SelfEnergyComplDMatrix;
 
  denomMatrComplInv.invert();

  Matrix< complex <double> > currentTMatr=_OmnesComplDMatrix*denomMatrComplInv*(*_Pvector);

  for (int i=0; i<currentTMatr.NumRows(); ++i){
    this->operator()(i,0)=currentTMatr(i,0);
  }

}

void FVectorOmnes::evalMatrix(const complex<double> mass, Spin OrbMom){

  _Kmatrix->evalMatrix(mass, OrbMom);
  _Pvector->evalMatrix(mass, OrbMom);


    for (int i=0; i<NumRows(); ++i){
    for (int j=0; j<NumRows(); ++j){
      _OmnesComplDMatrix(i,j) = _phpList->omnesMatrix()[i][j]->eval(mass*mass);
      _SelfEnergyComplDMatrix(i,j) = _phpList->selfEnergyMatrix()[i][j]->eval(mass*mass);
    }
  }
  Matrix< complex< double > > denomMatrComplInv = _idMatrix+(*_Kmatrix)*_SelfEnergyComplDMatrix;

  denomMatrComplInv.invert();

  Matrix< complex <double> > currentTMatr=_OmnesComplDMatrix*denomMatrComplInv*(*_Pvector);

  for (int i=0; i<currentTMatr.NumRows(); ++i){
    this->operator()(i,0)=currentTMatr(i,0);
  }

}

complex<double> FVectorOmnes::evalProjMatrix(const double mass, int index, Spin OrbMom){
  _Kmatrix->evalMatrix(mass, OrbMom);
  _Pvector->evalMatrix(mass, OrbMom);


      for (int i=0; i<NumRows(); ++i){
    for (int j=0; j<NumRows(); ++j){
      _OmnesComplDMatrix(i,j) = _phpList->omnesMatrix()[i][j]->eval(mass*mass);
      _SelfEnergyComplDMatrix(i,j) = _phpList->selfEnergyMatrix()[i][j]->eval(mass*mass);
    }
  }

  Matrix< complex< double > > denomMatrComplInv = _idMatrix+(*_Kmatrix)*_SelfEnergyComplDMatrix;
  denomMatrComplInv.invert();
  
  complex <double> result(0.,0.);
  for(int m=0;m<NumRows(); ++m){
    for(int n=0;n<NumRows(); ++n){
      result+=_OmnesComplDMatrix(index,m)*denomMatrComplInv(m,n)*(*_Pvector)(n,0);
    }
  }

  if(result.real() != result.real()){
    Alert << "result:\t" << result << endmsg;
    Alert << "mass:\t" << mass << "\tindex: " << index <<"\tOrbMom: " << OrbMom << endmsg;
    for(int i=0;i<NumRows(); ++i){
      Alert <<"(*_Pvector)(" << i << ",0): " << (*_Pvector)(i,0) << endmsg;
      Alert <<"denomMatrInv(index," << i << "): " << denomMatrComplInv(index,i) << endmsg;
    }
    exit(0);
  }

  return result; 
}



