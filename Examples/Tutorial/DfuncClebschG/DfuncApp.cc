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
#include "Examples/Tutorial/DfuncClebschG/DfuncHist.hh"

int main(int __argc,char *__argv[]){

  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){

    Info << "\nThis is a test application histograms the d functions depending on cos(theta) for a specific J lam1 and lam2\n"
	 << "\nd(J,lam1,lam2)\n"
         << "-2J   (default 4)\n"
         << "-2lam1   (default 0)\n"
         << "-2lam2   (default 0)"
         << endmsg;
    return 0;
  }

  std::string twoJStr="4";
  std::string twolam1Str="0";
  std::string twolam2Str="0";
  
  int optind = 1;
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-2J"){
      optind++;
      twoJStr = __argv[optind];
      found=true;
    }
    if (sw=="-2lam1"){
      optind++;
      twolam1Str = __argv[optind];
      found=true;
    }
    if (sw=="-2lam2"){
      optind++;
      twolam2Str = __argv[optind];
      found=true;
    }
    if (!found){
      Warning << "Unknown switch: " 
            << __argv[optind] << endmsg;
      optind++;
    }
    
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
    }

  std::stringstream twoJStrStr(twoJStr);
  int twoJ=0;
  twoJStrStr >> twoJ;

  std::stringstream twolam1StrStr(twolam1Str);
  int twolam1=0;
  twolam1StrStr >> twolam1;

  std::stringstream twolam2StrStr(twolam2Str);
  int twolam2=0;
  twolam2StrStr >> twolam2;

  Spin J( double(twoJ)/2. );
  Spin lam1( double(twolam1)/2. );
  Spin lam2( double(twolam2)/2. );

  Info <<"J= " << J << endmsg;
  Info <<"lam1= " << lam1 << endmsg;
  Info <<"lam2= " << lam2 << endmsg;

   DfuncHist dfuncHist(J, lam1, lam2);
  
 
  return 0;
}

