#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/Tutorial/LineShapes/TwoPoles.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);

  double mass1=1.35;
  double width1=0.3;

  double mass2=1.5;
  double width2=0.1; 

  
//   double mass1=1.2;
//   double width1=0.1;

//   double mass2=1.8;
//   double width2=0.1; 



  double delta=0.4;

  TwoPoles twoPoles(mass1, width1, mass2, width2, delta);

  return 0;
}

