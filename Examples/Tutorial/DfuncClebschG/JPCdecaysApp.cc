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

#include "Examples/Tutorial/DfuncClebschG/JPCdecays.hh"
#include "Examples/Tutorial/DfuncClebschG/DecAngularDist.hh"
#include "Examples/Tutorial/DfuncClebschG/SphericalHarmonicsPolts.hh"

#include "PwaUtils/DataUtils.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "ErrLogger/ErrLogger.hh"

bool findPCAndRemove(std::string& theString, int& theP, int& theC){

  bool result=true;

  std::string plusplusStr("++");
  std::string plusminusStr("+-");
  std::string minusplusStr("-+");
  std::string minusminusStr("--");
  std::string plus0Str("+0");
  std::string minus0Str("-0");

  size_t strfound=0;
  if ( (strfound=theString.find(plusplusStr)) != std::string::npos){
   theString.erase (strfound, strfound+2);
   theP=1;
   theC=1;
  }
  else  if ( (strfound=theString.find(plusminusStr)) != std::string::npos){
   theString.erase (strfound, strfound+2);
   theP=1;
   theC=-1;
  }
  else if ( (strfound=theString.find(minusplusStr)) != std::string::npos){
   theString.erase (strfound, strfound+2);
   theP=-1;
   theC=1;
  }
  else if ( (strfound=theString.find(minusminusStr)) != std::string::npos){
   theString.erase (strfound, strfound+2);
   theP=-1;
   theC=-1;
  }
  else if ( (strfound=theString.find(plus0Str)) != std::string::npos){
   theString.erase (strfound, strfound+2);
   theP=1;
   theC=0;
  }
  else if ( (strfound=theString.find(minus0Str)) != std::string::npos){
   theString.erase (strfound, strfound+2);
   theP=-1;
   theC=0;
  }
  else result=false;

  return result;
}



int main(int __argc,char *__argv[]){

  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){

    Info << "\nThis is a test application which calculates all allowed LS combinations for the decay a((2J)PC) -> b((2J)PC) + c((2J)PC)\n"
	 << "The switches are:\n\n"
         << "-mother  (default 2+-)\n\n" 
         << "-daughter1  (default 4--)\n\n"
         << "-daughter2   (default 2-+)"
         << endmsg;
    return 0;
  }

  std::string motherStr="2+-";
  std::string daughter1Str="4--";
  std::string daughter2Str="2-+";
  
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-mother"){
      optind++;
      motherStr = __argv[optind];
      found=true;
    }
    if (sw=="-daughter1"){
      optind++;
      daughter1Str = __argv[optind];
      found=true;
    }
    if (sw=="-daughter2"){
      optind++;
      daughter2Str = __argv[optind];
      found=true;
    }
    if (!found){
      Warning << "Unknown switch: " 
            << __argv[optind] << endmsg;
      optind++;
    }
    
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
    }

  Info << "the JPC combinations of the particles are:\n\n"
       << "mother:\t" << motherStr
       << "\ndaughter1:\t" << daughter1Str
       << "\ndaughter2:\t" << daughter2Str << endmsg;


  int motherP=0;
  int motherC=0;
  int daughter1P=0;
  int daughter1C=0;
  int daughter2P=0;
  int daughter2C=0;

  bool found=findPCAndRemove(motherStr, motherP, motherC);

  if (!found) {
    Alert << "JPCmother: " << motherStr << " not allowed!!!" << endmsg;
    exit(0);
  }

  std::stringstream mother2JStrStr(motherStr);
  double mother2J=0.;
  mother2JStrStr >> mother2J;
  Info << "mother:\tJ: " << mother2J << "/2\tP: " << motherP << "\tC: " << motherC << endmsg;

  found=findPCAndRemove(daughter1Str, daughter1P, daughter1C);
  if (!found) {
    Alert << "JPCdaughter1: " << daughter1Str << " not allowed!!!" << endmsg;
    exit(0);
  }
  std::stringstream daughter1_2JStrStr(daughter1Str);
  double daughter1_2J=0.;
  daughter1_2JStrStr >> daughter1_2J;
  Info << "daughter1:\tJ: " << daughter1_2J << "/2\tP: " << daughter1P << "\tC: " << daughter1C << endmsg;

  found=findPCAndRemove(daughter2Str, daughter2P, daughter2C);
  if (!found) {
    Alert << "JPCdaughter2: " << daughter2Str << " not allowed!!!" << endmsg;
    exit(0);
  }
  std::stringstream daughter2_2JStrStr(daughter2Str);
  double daughter2_2J=0.;
  daughter2_2JStrStr >> daughter2_2J;
  Info << "daughter2:\tJ: " << daughter2_2J << "/2\tP: " << daughter2P << "\tC: " << daughter2C << endmsg;


  
  std::shared_ptr<const jpcRes> motherJPC(new jpcRes(Spin(mother2J/2), motherP, motherC));
  std::shared_ptr<const jpcRes> daughter1JPC(new jpcRes(Spin(daughter1_2J/2), daughter1P, daughter1C)); 
  std::shared_ptr<const jpcRes> daughter2JPC(new jpcRes(Spin(daughter2_2J/2), daughter2P, daughter2C));
  motherJPC->print(std::cout);
  JPCdecays theDecays(motherJPC, daughter1JPC, daughter2JPC);
  theDecays.print(std::cout);

  DecAngularDist theDecAngularDist(motherJPC, daughter1JPC, daughter2JPC);

  SphericalHarmonicsPolts sphericalHarmonicsPolts(motherJPC->J+2); 
  return 0;
}

