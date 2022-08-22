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


PwaCovMatrix::PwaCovMatrix(std::map<std::string, std::map<std::string, double> >& theMaps){
  _n=theMaps.size();
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  std::map<std::string, double>::iterator it2;
  for(it1=theMaps.begin(); it1!=theMaps.end(); ++it1){
    for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      _covMatrix[it1->first][it2->first]=it2->second;
      if(it2->first != it1->first) _covMatrix[it2->first][it1->first]=it2->second;
    }
  }
  DiagonalIsValid();
  CheckCorrelationCoefficients();
}


PwaCovMatrix::PwaCovMatrix(ROOT::Minuit2::MnUserCovariance &theMinuitCovMatrix, 
			   const ROOT::Minuit2::MnUserParameters &theMinuitParameters)
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
	 InfoMsg << "Diagonal sqrt(cov) element (" << i << ", " << i << ") / param error: "
		 << sqrt(covValue) << " / " << errValue << endmsg;
       }
       else{
	 double correlationCoeff = covValue / sqrt( theMinuitCovMatrix(iCov, iCov) * theMinuitCovMatrix(jCov, jCov));
	 
	 if(fabs(correlationCoeff) > 0.99){
	   WarningMsg << "Correlation between parameter " << i << " and " << j << " = " << correlationCoeff << endmsg;
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
   if(it1 == _covMatrix.end()) {
     //     WarningMsg << "first element with name " << parameter1 << " not found!!!!" <<endmsg;
     return 0.;
   }

   it2 = (*it1).second.find(parameter2);
   if(it2 == (*it1).second.end()) {
     //     WarningMsg << "second element with name " << parameter2 << " not found!!!!" <<endmsg;
     return 0.;
   }

   return _covMatrix[parameter1][parameter2];
}

void PwaCovMatrix::SetElement(std::string parameter1, std::string parameter2, double val){

  std::map<std::string, std::map<std::string, double> >::iterator it1;
   std::map<std::string, double>::iterator it2;

   it1 = _covMatrix.find(parameter1);
   if(it1 == _covMatrix.end()) {
     WarningMsg << "first element with name " << parameter1 << " not found!!!!" <<endmsg;
     return;
     //     exit(1);
   }

   it2 = (*it1).second.find(parameter2);
   if(it2 == (*it1).second.end()) {
      WarningMsg << "second element with name " << parameter1 << " not found!!!!" <<endmsg;
     return;
     //     exit(1);
   }

   _covMatrix[parameter1][parameter2]=val;
   _covMatrix[parameter2][parameter1]=val;
}

void PwaCovMatrix::ResetAllElements(){
std::map<std::string, std::map<std::string, double> >::iterator it1;
  std::map<std::string, double>::iterator it2;

  for(it1=_covMatrix.begin(); it1!=_covMatrix.end(); ++it1){
    for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
    SetElement(it1->first, it2->first, 0.);
    }
  }
}

bool PwaCovMatrix::DiagonalIsValid(const ROOT::Minuit2::MnUserCovariance &theMinuitCovMatrix){

   bool result = true;
   for(unsigned int i=0; i<theMinuitCovMatrix.Nrow(); i++){
      if(theMinuitCovMatrix(i, i) < 0){
	 WarningMsg << "Covariance element (" << i << "," << i << ")"
		 << " = " << theMinuitCovMatrix(i,i) << " < 0" << endmsg;
	 result = false;
      }
   }
   return result;
}

bool PwaCovMatrix::DiagonalIsValid(){

   bool result = true;

   std::map<std::string, std::map<std::string, double> >::iterator it1;
   std::map<std::string, double>::iterator it2;
   
   for(it1=_covMatrix.begin(); it1!=_covMatrix.end(); ++it1){
     for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
       if(it1->first == it2->first){
	 if(it2->second < 0){
	   WarningMsg << "Covariance element (" << it1->first << " " << it2->first << ")"
		      << " = " << it2->second << " < 0" << endmsg;
	   result = false;
	 }
	 else{
	   InfoMsg << "Covariance element (" << it1->first << " " << it2->first << ")"
		   << " = " << it2->second << endmsg;
	 }
       }
     }
   }
   return result;
}

bool PwaCovMatrix::CheckCorrelationCoefficients(){
  bool result = true;

  std::map<std::string, std::map<std::string, double> >::iterator it1;
   std::map<std::string, double>::iterator it2;
   
   for(it1=_covMatrix.begin(); it1!=_covMatrix.end(); ++it1){
     for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
       if(it1->first != it2->first){
	 double denom=sqrt( GetElement(it1->first, it1->first) * GetElement(it2->first, it2->first));
	 if(std::abs(denom)>1.e-18){
	   double correlationCoeff = GetElement(it1->first, it2->first) / denom;
	   
	   if(fabs(correlationCoeff) > 0.999999){
	     WarningMsg << "Correlation between parameter " << it1->first << " and " << it2->first << " = " << correlationCoeff << " > 0.999999" << endmsg;
	     result = false;
	   }
	   else{
	      InfoMsg << "Correlation between parameter " << it1->first << " and " << it2->first << " = " << correlationCoeff << endmsg;
	   }
	 }
       }
     }
   }
   
   return result;
}

void PwaCovMatrix::printElements(){
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  std::map<std::string, double>::iterator it2;

  for(it1=_covMatrix.begin(); it1!=_covMatrix.end(); ++it1){
    for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      std::cout << it1->first << "\t" << it2->first 
		<< "\t" << it2->second << std::endl;
    }
  }
}
