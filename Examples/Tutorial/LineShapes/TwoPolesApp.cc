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

    Info << "\nThis is a test application for histogramming the mass shapes of two resonances obtained by the Breit-Wigner parameterization and by the T-matrix calculation\n" 
	 << "The switches are:\n\n"
         << "-mass1 (mass of the resonance1;default 1.318)\n\n" 
         << "-width1 (width1 of the resonance1;default 0.1)\n\n"
         << "-mass2 (mass of the resonance2;default 1.5)\n\n" 
         << "-delta (delta=mass range for histograms;default 0.4)\n\n" 
         << endmsg;
    return 0;
  }

  std::string mass1Str="1.318";
  std::string width1Str="0.1";
  std::string mass2Str="1.5";
  std::string width2Str="0.1";
  std::string deltaStr="0.1";

  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-mass1"){
      optind++;
      mass1Str = __argv[optind];
      found=true;
    }
    if (sw=="-width1"){
      optind++;
      width1Str = __argv[optind];
      found=true;
    }
    if (sw=="-mass2"){
      optind++;
      mass2Str = __argv[optind];
      found=true;
    }
    if (sw=="-width2"){
      optind++;
      width2Str = __argv[optind];
      found=true;
    }
    if (sw=="-delta"){
      optind++;
      deltaStr = __argv[optind];
      found=true;
    }
    if (!found){
      Warning << "Unknown switch: " 
	      << __argv[optind] << endmsg;
      optind++;
    }
  }


  std::stringstream mass1StrStr(mass1Str);
  double mass1=0.;
  mass1StrStr >> mass1;

  std::stringstream width1StrStr(width1Str);
  double width1=0.;
  width1StrStr >> width1;

  std::stringstream mass2StrStr(mass2Str);
  double mass2=0.;
  mass2StrStr >> mass2;

  std::stringstream width2StrStr(width2Str);
  double width2=0.;
  width2StrStr >> width2;

  std::stringstream deltaStrStr(deltaStr);
  double delta=0.;
  deltaStrStr >> delta;  

  TwoPoles twoPoles(mass1, width1, mass2, width2, delta);

  return 0;
}

