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

#include "qft++Extension/OmnesMatrixPiPiS3Channelwrhorho.hh"
#include "qft++Extension/UnityComplexFunction.hh"
#include "qft++Extension/ZeroComplexFunction.hh"

OmnesMatrixPiPiS3Channelwrhorho::OmnesMatrixPiPiS3Channelwrhorho() :
  AbsOmnesMatrix()
    ,_omnesMatrix(boost::extents[3][3])
  ,_selfEnergyMatrix(boost::extents[3][3])
{

  for(size_t i=0; i<3; ++i){
    for(size_t j=0; j<3; ++j){
      if(i==j){
        _omnesMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new UnityComplexFunction());
        _selfEnergyMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new UnityComplexFunction());
      }
      else{
        _omnesMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new ZeroComplexFunction());
        _selfEnergyMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new ZeroComplexFunction());
      }
    }
  }
}

OmnesMatrixPiPiS3Channelwrhorho::~OmnesMatrixPiPiS3Channelwrhorho(){
}

