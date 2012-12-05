
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>
#include "Examples/Tutorial/LineShapes/PiPiSWaveTMatrix.hh"
#include "Examples/Tutorial/LineShapes/PiPiSWavePVector.hh"
#include "Examples/Tutorial/LineShapes/KPiSWaveTMatrix.hh"

#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);

  PiPiSWaveTMatrix pipiSWaveTMatrix;
  PiPiSWavePVector piPiSWavePVector;
  KPiSWaveTMatrix kPiSWaveTMatrix;
  return 0;
}

