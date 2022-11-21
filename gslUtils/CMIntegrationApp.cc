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
    else if (sw=="-stepSizeImag"){
      optind++;
      stepSizeImag = atof(__argv[optind]);
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

  
  double mpole=0.90338253;
  double fpole=0.59405387;
  double mu=0.1349768;
  double m1=0.13957;
  double m2=0.13957;
  std::string oFileName("ChewMandelstaCalc.dat");
  std::ofstream oStream(oFileName.c_str());  

  InfoMsg << "...Starting with parameter..." << endmsg;
  InfoMsg << "mRealLow: " << mRealLow <<"\tmRealHigh: " << mRealHigh << "\tstepSizeReal: " << stepSizeReal << endmsg;
  InfoMsg << "mImagLow: " << mImagLow <<"\tmImagHigh: " << mImagHigh << "\tstepSizeImag: " << stepSizeImag << endmsg;

  
  CMIntegration cmIntegration(mpole, fpole, mu, m1, m2);

  for(double sqrtsReal=mRealLow; sqrtsReal<=mRealHigh; sqrtsReal+=stepSizeReal){
    for(double sqrtsImag=mImagHigh; sqrtsImag>=mImagLow; sqrtsImag-=stepSizeImag){
      std::complex<double> sqrts(sqrtsReal, sqrtsImag);
      std::complex<double> result=cmIntegration.integrate(sqrts*sqrts);
      oStream << sqrtsReal << "\t" << sqrtsImag << "\t" << result.real() << "\t" << result.imag() << std::endl;
      InfoMsg<<"CM(sqrt(s)=" << sqrts << ") = "<< result << endmsg;
    }
  }
  oStream.close(); 
  return 0;

}

