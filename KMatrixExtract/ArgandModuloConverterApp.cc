//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
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
#include <complex>
#include <map>
#include <getopt.h>
#include <iomanip>      // std::setprecision
#include <memory>

#include "KMatrixExtract/ArgandModuloConverter.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
  // ErrLogger::instance().setThreshold(logging::log_level::DEBUG);
  ErrLogger::instance().setThreshold(logging::log_level::INFO);

  if( __argc>1 && ( strcmp( __argv[1], "-h" ) == 0 ||
		strcmp( __argv[1], "--help" ) == 0 ) ) {
	InfoMsg << "USAGE:" << endmsg;
	InfoMsg << "--path: path to input root file file" << endmsg;
	InfoMsg << "--m1: mass of the first particle for the first php factor" << endmsg;
	InfoMsg << "--m2: mass of the second particle for the first php factor" << endmsg;
	InfoMsg << "--m3: mass of the third particle for the second php factor" << endmsg;
	InfoMsg << "--m4: mass of the fourth particle for the second php factor" << endmsg;
	return 0;
  }

  std::string pathToRootFile="";
  double mass1=0.135;
  double mass2=0.135; 
  double mass3=0.135;
  double mass4=0.135;

  for(;optind < (__argc-1); optind++){

	std::string ws = __argv[optind];

	if(ws[0]!='-' && ws[1]!='-'){
	  continue;
	}     
	else if (ws == "--path"){
	  pathToRootFile = __argv[optind+1];
	}
	else if (ws == "--m1"){
	  std::istringstream stream(__argv[optind+1]);
	  stream >> mass1;
	}
	else if (ws == "--m2"){
	  std::istringstream stream(__argv[optind+1]);
	  stream >> mass2;
	}
	else if (ws == "--m3"){
	  std::istringstream stream(__argv[optind+1]);
	  stream >> mass3;
	}
	else if (ws == "--m4"){
	  std::istringstream stream(__argv[optind+1]);
	  stream >> mass4;
	}
	else{
	  Alert << "Unknown switch: " << __argv[optind] << endmsg;
	  return 0;
	}
  }

  InfoMsg << "root file: " << pathToRootFile 
	  << "\nm1: " << mass1 << "\tm2: " << mass2 << "\tm3: " << mass3 << "\tm4: " << mass4
	  << endmsg;

  ArgandModuloConverter theConverter(pathToRootFile, mass1, mass2, mass3, mass4);

  theConverter.convertArgandToModulo(); 
  return 0;
}
