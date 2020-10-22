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

#include "PwaDynamics/PhaseSpaceList.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"
#include "qft++Extension/UnityComplexFunction.hh"
#include "qft++Extension/ZeroComplexFunction.hh"
#include "ErrLogger/ErrLogger.hh"

//PhaseSpaceList::PhaseSpaceList(vector<std::shared_ptr<AbsPhaseSpace> > phpVecs) :
//  _phpVecs(phpVecs)
//  ,_omnesMatrix(boost::extents[phpVecs.size()][phpVecs.size()])
//  ,_selfEnergyMatrix(boost::extents[phpVecs.size()][phpVecs.size()])
//{
//  for(size_t i=0; i<phpVecs.size(); ++i){
//    for(size_t j=0; j<phpVecs.size(); ++j){
//      if(i==j){
//       	_omnesMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new UnityComplexFunction());
//	_selfEnergyMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new UnityComplexFunction());
//      }
//      else{
//	_omnesMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new ZeroComplexFunction());
//	_selfEnergyMatrix[i][j] = std::shared_ptr<AbsComplexFunction>(new ZeroComplexFunction());
//      }
//    }
//  }
//}

PhaseSpaceList::PhaseSpaceList(vector<std::shared_ptr<AbsPhaseSpace> > phpVecs, boost::multi_array< std::shared_ptr<AbsComplexFunction> , 2> omnesMatr, boost::multi_array< std::shared_ptr<AbsComplexFunction> , 2> selfEnergyMatr) :
  _phpVecs(phpVecs)
  ,_omnesMatrix(omnesMatr)
  ,_selfEnergyMatrix(selfEnergyMatr)
{
  if(_omnesMatrix.num_elements() != phpVecs.size()*phpVecs.size()){
    Alert << "omnesMatrix number of elements: " << _omnesMatrix.num_elements() << " is not in agreement with phpSize*phpSize: " << phpVecs.size()*phpVecs.size() << endmsg;
    exit(1);
  }

    if(_selfEnergyMatrix.num_elements() != phpVecs.size()*phpVecs.size()){
    Alert << "omnesMatrix number of elements: " << _selfEnergyMatrix.num_elements() << " is not in agreement with phpSize*phpSize: " << phpVecs.size()*phpVecs.size() << endmsg;
    exit(1);
  }
}


PhaseSpaceList::~PhaseSpaceList(){

}
