//************************************************************************//
//									  //
//  Copyright 2022 Bertram Kopf (bertram@ep1.rub.de)		  //
//  	      	   - Ruhr-Universität Bochum 				  //
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
#include <fstream>
#include <sstream>
#include <vector>
#include <complex>
#include <getopt.h>

#include <memory>

#include "ErrLogger/ErrLogger.hh"


int main(int __argc,char *__argv[]){
  ErrLogger::instance().setThreshold(logging::log_level::INFO);
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){

    InfoMsg << "\nThis is an application for reading elements of the LUT files\n"
	    << "The switches are:\n\n"
	    << "-f (path to input LUT-file)\n\n" 
	    << endmsg;
    return 0;
  }

  std::cout << std::endl;
  std::cout << "___Starting..." << std::endl;

  std::string lutfilepath;
  
  // Read arguments and replace default values
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-f"){
      optind++;
      lutfilepath = __argv[optind];
      found=true;
    }
    // if (sw=="-m2"){
    //   optind++;
    //   p3Str = __argv[optind];
    //   found=true;
    // }
    if (!found){
      WarningMsg << "Unknown switch: " << __argv[optind] << endmsg;
      optind++;
    }
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
    }


  if (FILE *file = fopen(lutfilepath.c_str(), "r")) {
    fclose(file);
    InfoMsg << "input LUT file: " << lutfilepath.c_str() << endmsg;
  } else {
    Alert << "file: " << lutfilepath << " does not exist!!!!" << endmsg;
    exit(1);  
  }

  std::ifstream lutfile;
  lutfile.open(lutfilepath.c_str(), std::ios::binary | std::ifstream::in);
  //lutfile.open(lutfilepath.c_str(), std::ios::in | std::ios::binary);
  if (!lutfile.is_open()){
    Alert << "LUT: file not readable: " << lutfilepath << endmsg;
    exit(1);
  }

  //char ch;
  //in.read((char*)&ch, sizeof(ch));
  //std::cout << "Char: " << ch << std::endl;

	
  // get length of file:
  lutfile.seekg (0, std::ios::end);
  long size = lutfile.tellg();
  lutfile.seekg (0, std::ios::beg);

  double firstRe = 0.0;
  double lastRe = 0.0;
  double firstIm = 0.0;
  double lastIm = 0.0;
  int dSize = sizeof(lastRe);
  int m_dSize = dSize;
  long lines = size/(4*dSize);
  char *memblock = new char [dSize];
  lutfile.read(memblock, dSize);
  double *res=reinterpret_cast<double*>(memblock);

  firstRe = *res;
  lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  firstIm = *res;
  lutfile.seekg(size-4*dSize);
  lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  lastRe = *res;
  lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  lastIm = *res;

  std::complex<double> m_sLow(firstRe, firstIm);
  std::complex<double> m_sHigh(lastRe, lastIm);
  std::cout << "DEBUG-- sLow: " << m_sLow << " sHigh: " << m_sHigh << std::endl;

  lutfile.seekg(0, std::ios::beg);


  for(int i=0; i<lines; ++i){
    lutfile.read(memblock, dSize);
    double m_rel(*reinterpret_cast<double*>(memblock));
    //lutfile.seekg(i*4*dSize+dSize);
    lutfile.read(memblock, dSize);
    double m_imag(*reinterpret_cast<double*>(memblock));
    //lutfile.seekg(i*4*dSize+2*dSize);
    lutfile.read(memblock, dSize);
    double cm_rel(*reinterpret_cast<double*>(memblock));
    //lutfile.seekg(i*4*dSize+3*dSize);
    lutfile.read(memblock, dSize);
    double cm_imag(*reinterpret_cast<double*>(memblock));
    std::cout << m_rel << " " << m_imag << " " << cm_rel << " " << cm_imag << std::endl;
  }
  
  
  //int offset = 1;
  //int nReFoundInBlock = 1;
  //bool nReBlockFound = false;
  //  while((!nReBlockFound) && offset < size/(4*dSize)){
  //while((!nReBlockFound) && offset < size/(dSize)){
  //lutfile.seekg(offset*4*dSize);
  //  for (int i=0; i<4; ++i){ 
  //  lutfile.seekg(offset*i*dSize);
  //  lutfile.read(memblock, i*dSize);
  //  res=reinterpret_cast<double*>(memblock);
  //  std::cout << *res << " ";
  //  if(i==0) resloop=reinterpret_cast<double*>(memblock);
  //  if(abs(*res - firstRe) > 1.e-15) nReBlockFound = true;
  //  else nReFoundInBlock++;
  //  offset++;
  //  }
    
  //   std::cout << "\n" << std::endl;
    //if(abs(*resloop - firstRe) > 1.e-15) nReBlockFound = true;
    //else nReFoundInBlock++;
    //offset++;
  //  }
  //int m_nIm = nReFoundInBlock;
  //int  m_nRe = lines/m_nIm;
  
  // allocate memory:
  //char *  buffer = new char [length];
  // read data as a block:
  //lutfile.read(buffer,length);
  //lutfile.close();
  //std::cout << "Char: " << buffer << std::endl;

  //char ch;
  //lutfile.read((char*)&ch, sizeof(ch));
  //std::cout << "Char: " << ch << std::endl;
	
  return 0;

}

