#include <iostream>
#include <fstream>


#include "FitParamErrorMatrixStreamer.hh"
#include "ErrLogger/ErrLogger.hh"


FitParamErrorMatrixStreamer::FitParamErrorMatrixStreamer(  std::string errFile  ){
  std::ifstream inputStream;
  inputStream.open(errFile.c_str());
  if (!inputStream) {
   Alert << "can not open " << errFile ;
    exit(1);
  }
  while (!inputStream.eof()) {
    double err;
    inputStream >> err;
    if (!inputStream.fail()) {
      _matrixData.push_back(err);
     
    }
  }
  
}

