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

#include "Examples/Tutorial/LineShapes/VoigtShape.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]) {
  ErrLogger::instance().setThreshold(logging::log_level::DEBUG);
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ) {

    InfoMsg << "\nThis is a test application for histogramming the line shape of an Voigtian\n"
	    << "The switches are:\n\n"
	    << "-mass (mass of the resonance;default 1.318)\n\n" 
	    << "-width (width of the resonance;default 0.1)\n\n"
	    << "-sigma (resolution)\n\n"
	    << endmsg;
    return 0;
  }

  double mass0=0.782;
  double gamma=0.00849;
  double sigma=0.01;
  VoigtShape VoigtShape(mass0, gamma, sigma);

  return 0;
}

