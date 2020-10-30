//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
//          	   - Ruhr-Universit??t Bochum 				  //
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

#include "qft++Extension/OmnesMatrixPiPiS5Channel.hh"
#include "qft++Extension/UnityComplexFunction.hh"
#include "qft++Extension/ZeroComplexFunction.hh"
#include "qft++Extension/InputOmnesMatrixPiPiS5Channel.hh"
#include "qft++Extension/InputSelfEnergyMatrixPiPiS5Channel.hh"
#include "ErrLogger/ErrLogger.hh"
OmnesMatrixPiPiS5Channel::OmnesMatrixPiPiS5Channel() :
  AbsOmnesMatrix()
  ,_omnesMatrix(boost::extents[5][5])
  ,_selfEnergyMatrix(boost::extents[5][5])
{
  //  _omnesMatrix=boost::multi_array< std::shared_ptr<AbsComplexFunction> , 2> (boost::extents[5][5]);
  // _selfEnergyMatrix=boost::multi_array< std::shared_ptr<AbsComplexFunction> , 2> (boost::extents[5][5]);
  for(size_t i=0; i<5; ++i){
    for(size_t j=0; j<5; ++j){
      if(i<2 && j<2){
	_omnesMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new InputOmnesMatrixPiPiS5Channel(i,j));
        _selfEnergyMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new InputSelfEnergyMatrixPiPiS5Channel(i,j));
      }
      else if(i>=2 && j>=2 && i==j){
	_omnesMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new UnityComplexFunction());
	_selfEnergyMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new InputSelfEnergyMatrixPiPiS5Channel(i,j));
      }
      else{
        _omnesMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new ZeroComplexFunction());
        _selfEnergyMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new ZeroComplexFunction());
      }
    }
  }
}

OmnesMatrixPiPiS5Channel::~OmnesMatrixPiPiS5Channel(){
}

