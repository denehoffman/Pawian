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

#include "FitParams/PwaCovMatrix.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Minuit2/MnUserCovariance.h"
#include "Minuit2/MinuitParameter.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


PwaCovMatrix::PwaCovMatrix(){

}



PwaCovMatrix::PwaCovMatrix(ROOT::Minuit2::MnUserCovariance &theMinuitCovMatrix, 
			   const ROOT::Minuit2::MnUserParameters &theMinuitParameters,
			   fitParCol &theFitParams)
{
   DiagonalIsValid(theMinuitCovMatrix);

   _n = theMinuitCovMatrix.Nrow();
   unsigned int _nPar = theMinuitParameters.Params().size();
  
   int iCov=0;
   int jCov=0;

   for(unsigned int i=0; i<_nPar; i++){

      if(theMinuitParameters.Parameter(i).IsFixed())
	 continue;

      for(unsigned int j=0; j<_nPar; j++){

	 if(theMinuitParameters.Parameter(j).IsFixed())
	    continue;
 
	 double covValue = theMinuitCovMatrix(iCov, jCov);
	 double errValue = theMinuitParameters.Error(i);

	 if(iCov==jCov){
	    Info << "Diagonal sqrt(cov) element (" << i << ", " << i << ") / param error: "
		 << sqrt(covValue) << " / " << errValue << endmsg;
	 }
	 else{
	    double correlationCoeff = covValue / sqrt( theMinuitCovMatrix(iCov, iCov) * theMinuitCovMatrix(jCov, jCov));

	    if(fabs(correlationCoeff) > 0.99){
	       Warning << "Correlation between parameter " << i << " and " << j << " = " << correlationCoeff << endmsg;
	    }
	 }

	 _covMatrix[theMinuitParameters.GetName(i)][theMinuitParameters.GetName(j)] = covValue;

	 jCov++;
      }

      iCov++;
      jCov=0;
   }
}



double PwaCovMatrix::GetElement(std::string parameter1, std::string parameter2){

   std::map<std::string, std::map<std::string, double> >::iterator it1;
   std::map<std::string, double>::iterator it2;

   it1 = _covMatrix.find(parameter1);

   if(it1==_covMatrix.end()){
      return 0;
   }

   it2 = (*it1).second.find(parameter2);

   if(it2== (*it1).second.end()){
      return 0;
   }

   return _covMatrix[parameter1][parameter2];
}



bool PwaCovMatrix::DiagonalIsValid(const ROOT::Minuit2::MnUserCovariance &theMinuitCovMatrix){

   bool result = true;
   for(unsigned int i=0; i<theMinuitCovMatrix.Nrow(); i++){
      if(theMinuitCovMatrix(i, i) < 0){
	 Warning << "Covariance element (" << i << "," << i << ")"
		 << " = " << theMinuitCovMatrix(i,i) << " < 0" << endmsg;
	 result = false;
      }
   }
   return result;
}
