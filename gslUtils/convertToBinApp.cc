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

#include "ErrLogger/ErrLogger.hh"


int main(int __argc,char *__argv[]){
  ErrLogger::instance().setThreshold(logging::log_level::INFO);
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){

    InfoMsg << "\nThis is a application for calculating values of Chew-Mandelstam functions\n"
	    << "The switches are:\n\n"
	    << "-iFile (path to input file; default: ChewMandelstaCalc.dat)\n\n" 
	    << "-oFile (path to output file ; default: ChewMandelstaCalc.lut)\n\n"
	    << "-scalingFactor (scaling factor ; default: 1.0)\n\n" 
	    << endmsg;
    return 0;
  }

  // Set default values
  double scalingFactor=0.2;
  std::string iFile="./ChewMandelstaCalc.dat";
  std::string oFile="./ChewMandelstaCalc.lut";

  // Read arguments and replace default values
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    std::string sw = __argv[optind];
    if (sw=="-iFile"){
      optind++;
      iFile = __argv[optind];
    }
    else if (sw=="-oFile"){
      optind++;
      oFile = __argv[optind];
    }
    else if (sw=="-scalingFactor"){
      optind++;
      scalingFactor = atof(__argv[optind]);
    }
    else{
      WarningMsg << "Unknown switch: " << __argv[optind] << endmsg;
      optind++;
     }
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
  }

  //std::string oFileName("ChewMandelstaCalc.dat");
  //std::ofstream oStream(oFileName.c_str());  

  InfoMsg << "...Starting with conversion..." << endmsg;
  InfoMsg << "input file: " << iFile << endmsg;
  InfoMsg << "output file: " << oFile << endmsg;
  InfoMsg << "scaling factor: " << scalingFactor << endmsg;

  std::ifstream iStream(iFile.c_str());
  if (!iStream) {
    Alert << "can not open " << iFile << endmsg;
    exit(1);
  }

  std::ofstream oStream(oFile.c_str(), std::ios::out | std::ios::binary);
  //if (!oStream) {
  //   Alert << "can not open " << iFile << endmsg;
  //  exit(1);
  // }

  int counter=0;
  while (!iStream.eof()) {
    double mReal, mImag, cmReal, cmRealErr, cmImag, cmImagErr;
    iStream >> mReal >> mImag >> cmReal >> cmImag >> cmRealErr >> cmImagErr;
    //oStream << mReal << mImag << cmReal << cmImag << std::endl;
    //oStream << mReal << mImag << cmReal << cmImag;
    //if (counter>0) oStream << "\n";

    //oStream << mReal << mImag << cmReal << cmImag << "\n";
    cmReal*=scalingFactor;
    cmImag*=scalingFactor;
    
    oStream.write(reinterpret_cast<char*>(&mReal),sizeof(mReal));
    oStream.write(reinterpret_cast<char*>(&mImag),sizeof(mImag));
    oStream.write(reinterpret_cast<char*>(&cmReal),sizeof(cmReal));
    oStream.write(reinterpret_cast<char*>(&cmImag),sizeof(cmImag));
    counter++;
    std::cout << "mreal " << mReal << "  mimag " << mImag << "  counter: " << counter << std::endl;
  }

  oStream.close();
  return 0;

}

