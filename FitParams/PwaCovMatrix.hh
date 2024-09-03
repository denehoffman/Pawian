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

// PwaCovMatrix class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include <map>
#include <vector>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include "Minuit2/MnUserCovariance.h"
#include "Minuit2/MnUserParameters.h"


class PwaCovMatrix
{
  private:
   friend class boost::serialization::access;
   unsigned short _n;
   std::map<std::string, std::map<std::string, double> > _covMatrix;
   std::vector<double> _dataFromMnCovMatrix;

  public:
   PwaCovMatrix();
  PwaCovMatrix(std::map<std::string, std::map<std::string, double> >& theMaps);
   PwaCovMatrix(ROOT::Minuit2::MnUserCovariance &theMinuitCovMatrix, 
		const ROOT::Minuit2::MnUserParameters &theMinuitParameters);
   double GetElement(std::string parameter1, std::string parameter2);
   void SetElement(std::string parameter1, std::string parameter2, double val);
  void ResetAllElements();

   static bool DiagonalIsValid(const ROOT::Minuit2::MnUserCovariance &theMinuitCovMatrix);
  bool DiagonalIsValid();
  bool CheckCorrelationCoefficients();
  const std::vector<double>& dataFromMnCovMatrix() const {return _dataFromMnCovMatrix;}
  const unsigned int nRow() const {return _n;}

  void printElements();

   template<class Archive>
   void serialize(Archive & ar, const unsigned int version){
     ar & _n;
     ar & _covMatrix;
     if(version > 0){
       ar & _dataFromMnCovMatrix;
     }
   }

};

BOOST_CLASS_VERSION(PwaCovMatrix, 1)
