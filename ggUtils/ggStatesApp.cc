//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/AbsStates.hh"
#include "ggUtils/ggStates.hh"
#include "PwaUtils/DataUtils.hh"

int main(int __argc,char *__argv[]){

  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
		    || strcmp( __argv[1], "--help" ) == 0 ) ){

    std::cout << "test application which calculates all initial gamma gamma / GG states\n"
	      << "The application can be started with two flags: \n"
	      << "with -msg <errorLogMode>: choose the mode for the error logger \n"
	      << "with -jmax <value>: choose the maximum Spin J \n"
              << "e.g. type: ./ggStatesApp -jmax 6 -msg debugging \n"
	      << std::endl;
    return 0;
  }

  int optind=1;
  std::string msgModeStr="default";
  std::string jmaxStr="0";
  // decode arguments

  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-msg"){
      optind++;
      msgModeStr = __argv[optind];
      found=true;
    }
    if(sw=="-jmax"){
      optind++;
      jmaxStr = __argv[optind];
      found=true;
    }
    if (!found){
      WarningMsg << "Unknown switch: " 
	      << __argv[optind] << endmsg;
      optind++;
    }
    
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
  }

  if(msgModeStr == "debugging") ErrLogger::instance().setThreshold(logging::log_level::DEBUG);
  else if(msgModeStr == "trace") ErrLogger::instance().setThreshold(logging::log_level::INFO);
  else if(msgModeStr == "routine") ErrLogger::instance().setThreshold(logging::log_level::INFO);
  else if(msgModeStr == "warning") ErrLogger::instance().setThreshold(logging::log_level::WARN);
  else if(msgModeStr == "error")   ErrLogger::instance().setThreshold(logging::log_level::ERROR);
  else {
    ErrLogger::instance().setThreshold(logging::log_level::INFO);
    WarningMsg << "ErrorLogger not set -> Use mode 'INFO' " << endmsg;  
  }


 std::stringstream jmaxStrStr(jmaxStr);
 int jmax=0;
 jmaxStrStr >> jmax;

  InfoMsg << "jmax: " << jmax << endmsg;

  ggStates theggState(jmax);

  theggState.print(std::cout);

  return 0;
}

