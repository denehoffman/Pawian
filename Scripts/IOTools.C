//************************************************************************//
//                                                                        //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 - Ruhr-Universität Bochum                              //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

#pragma once

#include <cmath>
#include <iostream>

#include <fstream>
#include <sstream>
#include <string>

class IOTools
{
 public:
  static void addSysErrorToScatData(std::string iFileName, std::string oFileName, double error, double scaleVal1=1.){
    std::ifstream inputStream(iFileName.c_str());
    if(!inputStream.is_open()){
      std::cout << "Could not open input file: " << iFileName << std::endl;
      exit(0);
    }

    std::ofstream oStream(oFileName.c_str());

  while (!inputStream.eof()){
    double value1, value2, oldError;
    inputStream >> value1;
    inputStream >> value2;
    inputStream >> oldError;
   
    double newError=oldError+error;
    oStream << value1*scaleVal1 << "\t" << value2 << "\t" << newError << std::endl;     
  }

  inputStream.close();
  oStream.close();  
  }

  static void dumpFitparNameList(std::string iFileName, std::string oFileName, bool withMnParFix=true){
    std::string mnParFixStr="mnParFix = ";
    std::ifstream inputStream(iFileName.c_str());
    if(!inputStream.is_open()){
      std::cout << "Could not open input file: " << iFileName << std::endl;
      exit(0);
    }

    std::ofstream oStream(oFileName.c_str());

  while (!inputStream.eof()){
    std::string currentLine;
    getline (inputStream, currentLine);
    std::size_t pos = currentLine.find("\t"); 
    std::string parname = currentLine.substr(0,pos);
    if(parname=="") break;
    //    inputStream >> parname;
    std::cout << "parname: " << parname << std::endl;
    if(withMnParFix) oStream << mnParFixStr;  
    oStream << parname << std::endl;
  }

  inputStream.close();
  oStream.close();  
  }
};
