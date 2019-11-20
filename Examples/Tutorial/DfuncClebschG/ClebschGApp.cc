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


#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){

  ErrLogger::instance().setThreshold(logging::log_level::DEBUG);
  
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ) {

    InfoMsg << "\nThis is a test application calculates the Clebsch-Gordan coefficients "
	    << "C(j1, lambda_j1, j2, lambda_j2, J, lambda_J) with lambda_ji= -ji, -ji+1, ... , ji\n"
	    << "The switches are:\n\n"
	    << "-2j1  (default 1)\n\n" 
	    << "-2j2  (default 3)\n\n"
	    << "-2J   (default 4)"
	    << endmsg;
    return 0;
  }

  std::string twoj1Str="1";
  std::string twoj2Str="3";
  std::string twoJStr="4";
  
  int optind=1;
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-2j1"){
      optind++;
      twoj1Str = __argv[optind];
      found=true;
    }
    if (sw=="-2j2"){
      optind++;
      twoj2Str = __argv[optind];
      found=true;
    }
    if (sw=="-2J"){
      optind++;
      twoJStr = __argv[optind];
      found=true;
    }
    if (!found){
      WarningMsg << "Unknown switch: " << __argv[optind] << endmsg;
      optind++;
    }
    
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
    }

  std::stringstream twoj1StrStr(twoj1Str);
  unsigned int twoj1=0;
  twoj1StrStr >> twoj1;

  std::stringstream twoj2StrStr(twoj2Str);
  unsigned int twoj2=0;
  twoj2StrStr >> twoj2;

  std::stringstream twoJStrStr(twoJStr);
  unsigned int twoJ=0;
  twoJStrStr >> twoJ;

  Spin j1( double(twoj1)/2. );
  Spin j2( double(twoj2)/2. );
  Spin J( double(twoJ)/2. );

  if (J> j1+j2) {
    Alert << "J< j1+j2 !!!\t" << J << " < " << j1+j2 << endmsg;
    exit(0);
  }
  
  if (J< fabs(j1-j2)) {
    Alert << "J< |j1-j2| \t" << J << " < " << fabs(j1-j2) << endmsg;
    exit(0);
  }
  if ( (twoj1+twoj2)%2 != twoJ%2) {
    Alert << "j1= "<< j1 << " j2= " << j2 << "  cannot couple to J= " << J << endmsg;
    exit(0);
  }


  for (Spin lamj1=-j1;lamj1<=j1; ++lamj1) {
    for (Spin lamj2=-j2;lamj2<=j2; ++lamj2) {
      Spin lamJ=lamj1+lamj2;
      if (fabs(lamJ) > J) continue;
      InfoMsg << "C(j1=" << j1 << ", lamj1=" << lamj1 
	      << ", j2=" << j2 << ", lamj2=" << lamj2
	      << ", J=" << J << ", lamJ=" << lamJ << ")\t"
	      << Clebsch(j1, lamj1, j2, lamj2, J, lamJ)
	      << "\t  -> sqr= " 
	      << Clebsch(j1, lamj1, j2, lamj2, J, lamJ)*Clebsch(j1, lamj1, j2, lamj2, J, lamJ)
	      << endmsg;
    }
  }


  return 0;
}

