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

#include "Examples/Tutorial/LineShapes/TwoPoles.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);

  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){

    InfoMsg << "\nThis is a test application for histogramming the mass shapes of two resonances obtained by the Breit-Wigner parameterization and by the T-matrix calculation\n" 
	 << "The switches are:\n\n"
         << "-mass1 (mass of the resonance1;default 1.318)\n\n" 
         << "-width1 (width1 of the resonance1;default 0.1)\n\n"
         << "-mass2 (mass of the resonance2;default 1.5)\n\n" 
         << "-delta (delta=mass range for histograms;default 0.4)\n\n" 
         << endmsg;
    return 0;
  }

  double mass1=1.318;
  double width1=0.1;
  double mass2=1.5;
  double width2=0.1;
  double delta = 0.1;

  for(;optind < (__argc-1); optind++){

     std::string ws = __argv[optind];

     if(ws[0]!='-'){
       continue;
     }
     else if (ws == "-mass1" || ws == "--mass1"){
        std::istringstream mass1IStr(__argv[optind+1]);
	mass1IStr >> mass1;
     }
     else if (ws == "-width1" || ws == "--width1"){
        std::istringstream width1IStr(__argv[optind+1]);
	width1IStr >> width1;
     }
     else if (ws == "-mass2" || ws == "--mass2"){
        std::istringstream mass2IStr(__argv[optind+1]);
	mass2IStr >> mass2;
     }
     else if (ws == "-width2" || ws == "--width2"){
        std::istringstream width2IStr(__argv[optind+1]);
	width2IStr >> width2;
     }
     else if (ws == "-delta" || ws == "--delta"){
        std::istringstream deltaIStr(__argv[optind+1]);
	deltaIStr >> delta;
     }
     else{
        Alert << "Unknown switch: " << __argv[optind] << endmsg;
        return 0;
     }
  }

  TwoPoles twoPoles(mass1, width1, mass2, width2, delta);

  return 0;
}

