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

#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <memory>
#include <getopt.h>

#include "Examples/Tutorial/LineShapes/LUTPlot.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
  ErrLogger::instance().setThreshold(logging::log_level::DEBUG);
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ) {

    InfoMsg << "This is a test application for plotting various features of a "
	    << "BBUnstable phasespace factor lookup table\n"
	    << "The switches are:\n\n" 
	    << "-LUTFile path to the binary LUT file (string)\n\n"
	    << "-realGran granularity of the scan for the real part of sqrt(s) (double/GeV)\n\n"
	    << "-imagGran granularity of the scan for the imaginary part of sqrt(s) (double/GeV)\n\n"
	    << endmsg;

    return 0;
  }

  double realGranularity=0.01;
  double imagGranularity=0.01;
  std::string LUTFilePath = "";
  
  for(;optind < (__argc-1); optind++) {

    std::string ws = __argv[optind];
     if(ws[0]!='-'){
       continue;
     } 
    if (ws=="-LUTFile") {
      std::istringstream mass1IStr(__argv[optind+1]); 
      mass1IStr >> LUTFilePath; 
    }
    else if (ws=="-realGran") {
      std::istringstream mass2IStr(__argv[optind+1]); 
      mass2IStr >> realGranularity; 
    }
    else if (ws=="-imagGran") {
      std::istringstream mass3IStr(__argv[optind+1]); 
      mass3IStr >> imagGranularity; 
    }
    else{
      Alert << "Unknown switch: " << __argv[optind] << endmsg;
      return 1;
    }
  }
  if(LUTFilePath == "") {
    Alert << "You have to specify path to LUTFile with the -LUTFile switch! Exiting" << endmsg;
    return 1;
  }
  LUTPlot theLUTPlot(LUTFilePath, realGranularity, imagGranularity);  
  return 0;
}

