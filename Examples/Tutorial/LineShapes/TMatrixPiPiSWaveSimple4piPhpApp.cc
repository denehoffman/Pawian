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

#include "Examples/Tutorial/LineShapes/TMatrixPiPiSWaveSimple4piPhp.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);

  if( __argc>1 && ( strcmp( __argv[1], "-h" ) == 0 ||
                    strcmp( __argv[1], "--help" ) == 0 ) ){
    Info << "USAGE:" << endmsg;
    Info << "-s, --steps: number of steps in either direction in the complex energy plane" << endmsg;
    Info << "--maxImagMass: max imaginary part of the mass" << endmsg;
    Info << "--maxRealMass: max real part of the mass" << endmsg;
    Info << "--minImagMass: min imaginary part of the mass" << endmsg;
    Info << "--minRealMass: min reak part of the mass" << endmsg;
    return 0;
  }

  int numStepsForSheetScan = 500;

  std::vector<double> energyPlaneBorders;
  energyPlaneBorders.resize(4);
  energyPlaneBorders[0] = 0;    // Re min
  energyPlaneBorders[1] = -0.2; // Im min
  energyPlaneBorders[2] = 0;    // Re max
  energyPlaneBorders[3] = 0;    // Im max

  for(;optind < (__argc-1); optind++){

     std::string ws = __argv[optind];

     if(ws[0]!='-'){
       continue;
     }
      else if(ws == "--steps" || ws == "-s"){
	std::istringstream stream(__argv[optind+1]);
	stream >> numStepsForSheetScan;
     }
     else if(ws == "--minImagMass"){
	std::istringstream stream(__argv[optind+1]);
	stream >> energyPlaneBorders[1];
     }
     else if(ws == "--maxImagMass"){
	std::istringstream stream(__argv[optind+1]);
	stream >> energyPlaneBorders[3];
     }
     else if(ws == "--minRealMass"){
	std::istringstream stream(__argv[optind+1]);
	stream >> energyPlaneBorders[0];
     }
     else if(ws == "--maxRealMass"){
	std::istringstream stream(__argv[optind+1]);
	stream >> energyPlaneBorders[2];
     }
     else{
	Alert << "Unknown switch: " << __argv[optind] << endmsg;
	return 0;
     }
  }

  MatrixPiPiSWaveSimple4piPhp tMatrixPiPiSWaveSimple4piPhp(numStepsForSheetScan, energyPlaneBorders);

  return 0;
}

