//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#include "qft++Extension/OmnesMatrixFactory.hh"
#include "qft++Extension/AbsOmnesMatrix.hh"
#include "qft++Extension/OmnesMatrixGeneral.hh"
#include "ConfigParser/globalParser.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
  ErrLogger::instance().setThreshold(logging::log_level::DEBUG);

  if( __argc>1 && ( strcmp( __argv[1], "-h" ) == 0 ||
                     strcmp( __argv[1], "--help" ) == 0 ) ){
     InfoMsg << "USAGE:" << endmsg;
     InfoMsg << "-OmnesName: name/key of OmnesMatrix" << endmsg;
     return 0;
   }

  std::string omnesKey("");
  
  for(;optind < (__argc-1); optind++) {

    std::string ws = __argv[optind];
     if(ws[0]!='-') {
       continue;
     } 
    if (ws=="-OmnesName") {
      std::istringstream omnesKeyStr(__argv[optind+1]); 
      omnesKeyStr >> omnesKey; 
    }
    else {
      Alert << "Unknown switch: " << __argv[optind] << endmsg;
      return 0;
    }
  }
  std::shared_ptr<AbsOmnesMatrix> omnesMatrPtr=OmnesMatrixFactory::instance()->getOmnesMatrixPointer(omnesKey);

  OmnesMatrixGeneral omnesMatGeneral(omnesMatrPtr);
  
  return 0;
}

