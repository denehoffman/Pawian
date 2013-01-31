// PwaCovMatrix class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#include "PwaUtils/PwaCovMatrix.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Minuit2/MnUserCovariance.h"
#include "Minuit2/MinuitParameter.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


PwaCovMatrix::PwaCovMatrix(){

}



PwaCovMatrix::PwaCovMatrix(ROOT::Minuit2::MnUserCovariance &theMinuitCovMatrix, 
			   ROOT::Minuit2::MnUserParameters &theMinuitParameters,
			   fitParams &theFitParams)
{
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
      //Alert << "Covariance line not found: " << parameter1 << endmsg;
      return 0;
   }

   it2 = (*it1).second.find(parameter2);

   if(it2== (*it1).second.end()){
      //Alert << "Covariance column not found: " << parameter2 << endmsg;
      return 0;
   }

   return _covMatrix[parameter1][parameter2];
}

