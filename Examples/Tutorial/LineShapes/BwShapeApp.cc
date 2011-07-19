#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/Tutorial/LineShapes/BwShape.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){

    Info << "\nThis is a test application for histogramming the Argand plot, the phase shift and the mass shape of the relativistic Breit-Wigner function with Blatt-Weisskopf barrier factors\n"
	 << "The switches are:\n\n"
         << "-mass (mass of the resonance;default 1.318)\n\n" 
         << "-width (width of the resonance;default 0.1)\n\n"
         << "-mDec1 (mass of the first decay particle; default 0.13498\n\n"
         << "-mDec2 (mass of the second decay particle; default 0.547853\n\n"
         << "-Lmax (histogams of mass shapes; default 3\n\n"
         << "-DeltaM (range of histograms +-DeltaM; default 0.4\n\n"  
         << endmsg;
    return 0;
  }
  std::string massStr="1.318";
  std::string widthStr="0.1";
  std::string mDec1Str="0.13498";
  std::string mDec2Str="0.547853";
  std::string LmaxStr="3";
  std::string DeltaMStr="0.4";

  
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-mass"){
      optind++;
      massStr = __argv[optind];
      found=true;
    }
    if (sw=="-width"){
      optind++;
      widthStr = __argv[optind];
      found=true;
    }
    if (sw=="-mDec1"){
      optind++;
      mDec1Str = __argv[optind];
      found=true;
    }
    if (sw=="-mDec2"){
      optind++;
      mDec2Str = __argv[optind];
      found=true;
    }
    if (sw=="-Lmax"){
      optind++;
      LmaxStr = __argv[optind];
      found=true;
    }
    if (sw=="-DeltaM"){
      optind++;
      DeltaMStr = __argv[optind];
      found=true;
    }

    if (!found){
      Warning << "Unknown switch: " 
            << __argv[optind] << endmsg;
      optind++;
    }
    
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
    }


  std::stringstream massStrStr(massStr);
  double mass=0.;
  massStrStr >> mass;

  std::stringstream widthStrStr(widthStr);
  double width=0.;
  widthStrStr >> width;

  std::stringstream mDec1StrStr(mDec1Str);
  double mDec1=0.;
  mDec1StrStr >> mDec1;

  std::stringstream mDec2StrStr(mDec2Str);
  double mDec2=0.;
  mDec2StrStr >> mDec2;

  std::stringstream LmaxStrStr(LmaxStr);
  unsigned int Lmax=0;
  LmaxStrStr >> Lmax;

  std::stringstream DeltaMStrStr(DeltaMStr);
  double DeltaM=0;
  DeltaMStrStr >> DeltaM;

  BwShape bwShape(mass, width, mDec1, mDec2, Lmax, DeltaM);

  return 0;
}

