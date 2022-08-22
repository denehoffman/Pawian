//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
//          	   - Ruhr-Universit??t Bochum 				  //
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
#include <string>
#include <cstdlib>
#include <vector>
#include <map>
#include <tuple>
#include <iterator>
#include <memory>

#include "FitParams/PwaCovMatrix.hh"
#include "Utils/ErrLogUtils.hh"
#include "ErrLogger/ErrLogger.hh"


int main(int __argc,char *__argv[]){
  ErrLogger::instance().setThreshold(logging::log_level::DEBUG);

  if( __argc>1 && ( strcmp( __argv[1], "-h" ) == 0 ||
                    strcmp( __argv[1], "--help" ) == 0 ) ){
    InfoMsg << "USAGE:" << endmsg;
    InfoMsg << "-f, --file: path to serialization file" << endmsg;
    exit(0);
  }

  std::string pathToSerializationFile = "./serializedOutput.dat";
  for(;optind < (__argc-1); optind++){

     std::string ws = __argv[optind];if(ws[0]!='-'){
       continue;
     }
     else if (ws == "--file" || ws == "-f"){
        pathToSerializationFile = __argv[optind+1];
     }
 }
  InfoMsg << "path to pathToSerializationFile: " << pathToSerializationFile << endmsg;

  std::ifstream serializationStream(pathToSerializationFile.c_str());

  if(!serializationStream.is_open()){
        Alert << "Could not open serialization file." << endmsg;
        exit(0);
  }

  std::shared_ptr<PwaCovMatrix> thePwaCovMatrix = std::shared_ptr<PwaCovMatrix>(new PwaCovMatrix);
  boost::archive::text_iarchive boostInputArchive(serializationStream);

  boostInputArchive >> *thePwaCovMatrix;
  thePwaCovMatrix->printElements();

  //bool isDiagonalValid=thePwaCovMatrix->DiagonalIsValid();
  //InfoMsg << "isDiagonalValid: " << isDiagonalValid << endmsg;

  bool checkCorrelationCoefficients=thePwaCovMatrix->CheckCorrelationCoefficients();
  InfoMsg << "checkCorrelationCoefficients: " << checkCorrelationCoefficients << endmsg;
  return 1;
}


