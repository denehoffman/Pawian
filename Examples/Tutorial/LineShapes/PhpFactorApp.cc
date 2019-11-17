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

#include "Examples/Tutorial/LineShapes/PhpFactor.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
  ErrLogger::instance().setThreshold(logging::log_level::DEBUG);
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){

    InfoMsg << "This is a test application for scanning the real and imaginary part of the phase space factor\n"
	 << "The switches are:\n\n" 
         << "-m1 mass of the first decay particle (default: 0.135)\n\n"
         << "-m2 mass of the second decay particle (default: 0.547)\n\n"
         << "-mmax maximal mass of the scan (default: 2.0)\n"  
         << endmsg;

    return 0;
  }

  double mass1=0.135;
  double mass2=0.547;
  double massMax=2.0;
  
  for(;optind < (__argc-1); optind++) {

    std::string ws = __argv[optind];
     if(ws[0]!='-'){
       continue;
     } 
    if (ws=="-m1"){
      std::istringstream mass1IStr(__argv[optind+1]); 
      mass1IStr >> mass1; 
    }
    else if (ws=="-m2"){
      std::istringstream mass2IStr(__argv[optind+1]); 
      mass2IStr >> mass2;
    }
    else if (ws=="-mmax"){
      std::istringstream massMaxIStr(__argv[optind+1]); 
      massMaxIStr >> massMax;
    }
    else{
      Alert << "Unknown switch: " << __argv[optind] << endmsg;
      return 0;
    }
  }

  PhpFactor thePhpFactor(mass1, mass2, massMax);  
  return 0;
}

