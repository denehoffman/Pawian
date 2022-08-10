//************************************************************************//
//									  //
//  Copyright 2022 Bertram Kopf (bertram@ep1.rub.de)			  //
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
    InfoMsg << "-ifile: path to input serialization file" << endmsg;
    InfoMsg << "-ofile: path to output serialization file" << endmsg;
    InfoMsg << "-iascii: path to input ascii file with new parameters" << endmsg;
    InfoMsg << "-reset: reset remaining content" << endmsg;
    exit(0);
  }

  std::string pathToSerializationFile = "./serializedOutput.dat";
  std::string pathToNewSerializationFile = "./serializedNew.dat";
  std::string pathToAsciiFile = "./serializedNew.dat";
  
   while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-ifile") {
      optind++;
      pathToSerializationFile = __argv[optind];
      found=true;
    }
    if (sw=="-ofile") {
      optind++;
      pathToNewSerializationFile = __argv[optind];
      found=true;
    }
    if (sw=="-iascii") {
      optind++;
      pathToAsciiFile = __argv[optind];
      found=true;
    }
    if (!found) {
      Alert << "Unknown switch: " << __argv[optind] << endmsg;
      optind++;
      exit(1);
    }
    
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
    
   }
  
  InfoMsg << "path to pathToSerializationFile: " << pathToSerializationFile << endmsg;
  InfoMsg << "path to pathToNewSerializationFile: " << pathToNewSerializationFile << endmsg;
  InfoMsg << "path to pathToAsciiFile: " << pathToAsciiFile << endmsg;
  

  
  std::ifstream serializationStream(pathToSerializationFile.c_str());

  if(!serializationStream.is_open()){
        Alert << "Could not open serialization file." << endmsg;
        exit(0);
  }

  std::shared_ptr<PwaCovMatrix> thePwaCovMatrix = std::shared_ptr<PwaCovMatrix>(new PwaCovMatrix);
  boost::archive::text_iarchive boostInputArchive(serializationStream);
  boostInputArchive >> *thePwaCovMatrix;
  //thePwaCovMatrix->printElements();

  thePwaCovMatrix->ResetAllElements();
  
  std::ifstream asciiStream(pathToAsciiFile.c_str());

  if(!asciiStream.is_open()){
        Alert << "Could not open input ascii file." << endmsg;
        exit(0);
  }

  std::string firstName;
  std::string secondName;
  double value;
  while (!asciiStream.eof()){
    std::string currentLine;
    getline (asciiStream, currentLine);
    std::istringstream  lineiString(currentLine);
    lineiString >> firstName;
    lineiString	>> secondName;
    lineiString >> value;
    InfoMsg << firstName << "\t" << secondName << "\t" << value << endmsg;
    thePwaCovMatrix->SetElement(firstName, secondName, value);
    //    std::size_t pos = currentLine.find("\t"); 
    //    std::string currentname = currentLine.substr(0,pos);
    //    if(currentname=="") break;
    //     std::cout << "currentname: " << currentname << std::endl;
     //   if(withMnParFix) oStream << mnParFixStr;  
     //    oStream << parname << std::endl;
  }
  thePwaCovMatrix->printElements();

  // dump new serialization file
  std::ostringstream serializationFileName;
  std::ofstream serialization_oStream(pathToNewSerializationFile.c_str());
  boost::archive::text_oarchive boostOutputArchive(serialization_oStream);
  boostOutputArchive << *thePwaCovMatrix;
  return 1;
}


