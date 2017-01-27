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
#include <getopt.h>

#include <memory>

#include "Examples/Tutorial/LineShapes/FlatteShape.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){

    InfoMsg << "\nThis is a test application for histogramming the Flatte shape\n"
	 << "The switches are:\n\n"
         << "-ptype particle type (so far: a0, f0) (default a0)\n\n" 
         << "-g11 g(pi,eta) for a0; or g(pi,pi) for f0 (default 0.5857)\n\n"
         << "-g22 g(K,K) for a0 and f0 (default 1.86243)\n\n"
         << "-mass for a0 or f0  (default 0.982)\n"  
         << endmsg;
    return 0;
  }

  std::string ptypeStr="a0";
  std::string g11Str="0.5857";
  std::string g22Str="1.86243";
  std::string massStr="0.982";
  
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-ptype"){
      optind++;
      ptypeStr = __argv[optind];
      found=true;
    }
    if (sw=="-g11"){
      optind++;
      g11Str = __argv[optind];
      found=true;
    }
    if (sw=="-g22"){
      optind++;
      g22Str = __argv[optind];
      found=true;
    }
    if (sw=="-mass"){
      optind++;
      massStr = __argv[optind];
      found=true;
    }
    if (!found){
      WarningMsg << "Unknown switch: " 
            << __argv[optind] << endmsg;
      optind++;
    }
    
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
    }

  std::stringstream g11StrStr(g11Str);
  double g11=0.;
  g11StrStr >> g11;

  std::stringstream g22StrStr(g22Str);
  double g22=0.;
  g22StrStr >> g22;

  std::stringstream massStrStr(massStr);
  double mass=0.;
  massStrStr >> mass;
  

//   std::string resStr="a0";
//  std::string resStr="f0";


  FlatteShape flatteShape(ptypeStr, g11, g22);

  return 0;
}

