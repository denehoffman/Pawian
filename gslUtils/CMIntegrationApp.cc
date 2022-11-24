//************************************************************************//
//									  //
//  Copyright 2023 Bertram Kopf (bertram@ep1.rub.de)    		  //
//  Copyright 2023 Meike Kuessner (mkuessner@ep1.rub.de)                  //
//  	      	   - Ruhr-Universität Bochum 				  //
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

#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>

#include <vector>
#include <getopt.h>

#include <memory>
#include <complex>

#include "gslUtils/CMIntegration.hh"
#include "gslUtils/CMIntegrationDudek.hh"
#include "gslUtils/CMIntegrationReid.hh"
#include "ErrLogger/ErrLogger.hh"


int main(int __argc,char *__argv[]){
  ErrLogger::instance().setThreshold(logging::log_level::INFO);
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){

    InfoMsg << "\nThis is a application for calculating values of Chew-Mandelstam functions\n"
	    << "The switches are:\n\n"
	    << "-mRealLow (starting point of real mass; default 0.2)\n\n" 
	    << "-mRealHigh (terminal point of real mass ; default 3.0)\n\n"
            << "-stepSizeReal (step size real mass ; default 0.01)\n\n"
            << "-mImagLow (terminal point of imag mass; default -0.4)\n\n"
            << "-mRealHigh (starting point of imag mass ; default 0.0)\n\n"
      	    << "-stepSizeImag (step size imag mass ; default 0.01)\n\n"
	    << "-dumpErr (dump uncertaities for real and imagianry part of the CM; default false)\n\n"
	    << "-CMfunction (Dudek or Reid; default Dudek)\n\n"      
            << "-offsetMethod (offset method for Dudek inegration; default Reid)\n\n"
	    << "-mPole (default 0.9033435)\n\n"
	    << "-fPole (default 0.59399735)\n\n"
            << "-mu (default 0.13957)\n\n"
            << "-m1 (default 0.1349768)\n\n"
	    << "-m2 (default 0.13957)\n\n" 
	    << endmsg;
    return 0;
  }

  // Set default values
  double mRealLow=0.2;
  double mRealHigh=3.0;
  double stepSizeReal=0.01;
  double mImagLow=-0.4;
  double  mImagHigh=0.;
  double stepSizeImag=0.01;
  bool dumpWithErrors=false;
  std::string CMfunction="Dudek";
  std::string offsetMethod="Reid";
  double mpole=0.9033435;
  double fpole=0.59399735;
  double mu=0.13957;
  double m1=0.1349768;
  double m2=0.13957;

  // Read arguments and replace default values
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    std::string sw = __argv[optind];
    if (sw=="-mRealLow"){
      optind++;
      mRealLow = atof(__argv[optind]);
    }
    else if (sw=="-mRealHigh"){
      optind++;
      mRealHigh = atof(__argv[optind]);
    }
    else if (sw=="--stepSizeReal"){
      optind++;
      stepSizeReal = atof(__argv[optind]);
    }
    else if (sw=="-mImagLow"){
      optind++;
      mImagLow = atof(__argv[optind]);
    }
    else if (sw=="-mImagHigh"){
      optind++;
      mImagHigh = atof(__argv[optind]);
      //if(mImagHigh>0) mImagHigh=-mImagHigh;
    }
    else if (sw=="-dumpErr"){
      optind++;
      std::string currentStr = __argv[optind];
      dumpWithErrors = (currentStr == "1");
      if(!dumpWithErrors) dumpWithErrors = (currentStr == "true"); 
    }
    else if (sw=="-stepSizeImag"){
      optind++;
      stepSizeImag = atof(__argv[optind]);
    }
    else if (sw=="-offsetMethod"){
      optind++;
      offsetMethod = __argv[optind];
    }
    else if (sw=="-mPole"){
      optind++;
      mpole = atof(__argv[optind]);
    }
    else if (sw=="-fPole"){
      optind++;
      fpole = atof(__argv[optind]);
    }
    else if (sw=="-mu"){
      optind++;
      mu = atof(__argv[optind]);
    }
    else if (sw=="-m1"){
      optind++;
      m1 = atof(__argv[optind]);
    }
    else if (sw=="-m2"){
      optind++;
      m1 = atof(__argv[optind]);
    }
    else if (sw=="-CMfunction"){
      optind++;
      CMfunction = __argv[optind];
    }
    else{
      WarningMsg << "Unknown switch: " << __argv[optind] << endmsg;
      optind++;
     }
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
  }

  if(mImagHigh<mImagLow){
    Alert << "mImagHigh<mImagLow" << mImagHigh << " < " << mImagLow << endmsg;
    exit(1);
  }
  else if (mRealLow>mRealHigh){
    Alert << "mRealLow>mRealHigh" << mRealLow << " > " << mRealHigh << endmsg;
    exit(1);
  }

  
  // double mpole=0.9033435;
  // double fpole=0.59399735;
  // double mu=0.1349768;
  // double m1=0.13957;
  // double m2=0.13957;
  std::string oFileName("ChewMandelstaCalc.dat");
  std::ofstream oStream(oFileName.c_str());  

  InfoMsg << "...Starting with parameter..." << endmsg;
  InfoMsg << "mRealLow: " << mRealLow <<"\tmRealHigh: " << mRealHigh << "\tstepSizeReal: " << stepSizeReal << endmsg;
  InfoMsg << "mImagLow: " << mImagLow <<"\tmImagHigh: " << mImagHigh << "\tstepSizeImag: " << stepSizeImag << endmsg;
  InfoMsg << "dumpWithErrors: " << dumpWithErrors << endmsg;
  InfoMsg << "CMfunction: " << CMfunction << endmsg;
  InfoMsg << "offsetMethod: " << offsetMethod << endmsg;
  InfoMsg << "mPole: " << mpole << endmsg;
  InfoMsg << "fPole: " << fpole << endmsg;
  InfoMsg << "mu: " << mu << endmsg;
  InfoMsg << "m1: " << m1 << endmsg;
  InfoMsg << "m2: " << m2 << endmsg;
  
  CMunstable_params theParams(mpole, fpole, mu, m1, m2);  

  CMIntegration* cmIntegration=0;
  if (CMfunction=="Dudek") cmIntegration = new CMIntegrationDudek(mpole, fpole, mu, m1, m2, offsetMethod);
  else if (CMfunction=="Reid") cmIntegration = new CMIntegrationReid(mpole, fpole, mu, m1, m2);
  else{
    Alert << "CMfunction with the name " << CMfunction << " is not supported!!!" << endmsg;
    exit(0);
  }

  for(double sqrtsReal=mRealLow; sqrtsReal<=mRealHigh; sqrtsReal+=stepSizeReal){
    for(double sqrtsImag=mImagHigh; sqrtsImag>=mImagLow; sqrtsImag-=stepSizeImag){
      std::complex<double> currentResult;
      std::complex<double> currentResultErr;
      std::complex<double> sqrts(sqrtsReal, sqrtsImag);
      //      std::complex<double> result=cmIntegration.integrate(sqrts*sqrts, currentResult, currentResultErr);
      cmIntegration->integrate(sqrts*sqrts, currentResult, currentResultErr);
      if(dumpWithErrors) oStream << sqrtsReal << "\t" << sqrtsImag << "\t" << currentResult.real() << "\t" << currentResult.imag() << "\t"<<  currentResultErr.real()  << "\t"<<  currentResultErr.imag() << std::endl;
      else oStream << sqrtsReal << "\t" << sqrtsImag << "\t" << currentResult.real() << "\t" << currentResult.imag() << std::endl;
      InfoMsg<<"CM(sqrt(s)=" << sqrts << ") = "<< currentResult << endmsg;
    }
  }
  oStream.close();
  return 0;

}

