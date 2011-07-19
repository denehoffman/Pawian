#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/Tutorial/LineShapes/FlatteShape.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){

    Info << "\nThis is a test application for histogramming the Flatte shape\n"
	 << "The switches are:\n\n"
         << "-ptype particle type (so far: a0, f0) (default a0)\n\n" 
         << "-g11 g(pi,eta) for a0; or g(pi,pi) for f0 (default 0.5857)\n\n"
         << "-g22 g(K,K) for a0 and f0 (default 1.86243)\n\n"
         << "-mass for a0 or f0  (default 0.982)\n"  
         << endmsg;
    return 0;
  }

  std::string ptypeStr="a0";
  std::string g11Str="0.5857";
  std::string g22Str="1.86243";
  std::string massStr="0.982";
  
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-ptype"){
      optind++;
      ptypeStr = __argv[optind];
      found=true;
    }
    if (sw=="-g11"){
      optind++;
      g11Str = __argv[optind];
      found=true;
    }
    if (sw=="-g22"){
      optind++;
      g22Str = __argv[optind];
      found=true;
    }
    if (sw=="-mass"){
      optind++;
      massStr = __argv[optind];
      found=true;
    }
    if (!found){
      Warning << "Unknown switch: " 
            << __argv[optind] << endmsg;
      optind++;
    }
    
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
    }

  std::stringstream g11StrStr(g11Str);
  double g11=0.;
  g11StrStr >> g11;

  std::stringstream g22StrStr(g22Str);
  double g22=0.;
  g22StrStr >> g22;

  std::stringstream massStrStr(massStr);
  double mass=0.;
  massStrStr >> mass;
  

//   std::string resStr="a0";
//  std::string resStr="f0";


  FlatteShape flatteShape(ptypeStr, g11, g22);

  return 0;
}

