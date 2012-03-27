#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "Examples/EtacToapi0Fit/EtacToapi0EventList.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Hist.hh"
#include "Examples/EtacToapi0Fit/MEtacToapi0Fcn.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Data.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Lh.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

//#include "Minuit2/MnUserTransformation.h"
using namespace ROOT::Minuit2;

int main(int __argc,char *__argv[]){

  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
		    || strcmp( __argv[1], "--help" ) == 0 ) ){

    std::cout << "This application is a simple PWA fit for the decay chain\n"
	      << "eta_c -> intermediate + pi0  ; intermediate -> pi0 + eta\n"
	      << "It makes use of the Covariant Tensor Formalism\n"
	      << "The fit determines whether the intermediate resonance is a spin 0,1 or 2 particle\n"
	      << "In addition the mass and width of the intermediate resonance will be fitted with a simple Breit-Wigner\n"
              << "To start the application with data containing  an intermediate resonance with Spin=0 and mass 0.98 GeV, type: ./MEtacToapi0FitApp -d 0\n"
              << "To start the application with data containing  an intermediate resonance with Spin=2 and mass 1.32 GeV, type: ./MEtacToapi0FitApp -d 2\n"
              << "with the flag -msg <errorLogMode> you can choose the mode for the error logger\n"
	      << "i.e. with './MEtacToapi0FitApp -d 2 -msg debugging' you start the fit with Spin=2 in the debugging mode for the error logger\n"  
	      << std::endl;
    return 0;
  }

  
  int optind=1;
  std::string dataSpinStr="";
  std::string msgModeStr="default";
  // decode arguments
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-d") {
      optind++;
      dataSpinStr = __argv[optind];
      found=true;
    }
    if (sw=="-msg"){
      optind++;
      msgModeStr = __argv[optind];
      found=true;
    }
    if (!found){
      Warning << "Unknown switch: " 
            << __argv[optind] << endmsg;
      optind++;
    }
    
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
  }

  if(msgModeStr == "debugging") ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  else if(msgModeStr == "trace") ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
  else if(msgModeStr == "routine") ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
  else if(msgModeStr == "warning") ErrLogger::instance()->setLevel(log4cpp::Priority::WARN);
  else if(msgModeStr == "error")   ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR);
  else {
    ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
    Warning << "ErrorLogger not (properly) set -> Use mode 'DEBUG' " << endmsg;  
  }



 std::stringstream dataSpinStrStr(dataSpinStr);
 int dataSpin=2;
 dataSpinStrStr >> dataSpin ;

 Info << "dataSpin: " << dataSpin << endmsg;

 boost::shared_ptr<const EtacToapi0EventList> theEvtListPtr(new EtacToapi0EventList(dataSpin));
 
 boost::shared_ptr<EtacToapi0Lh> theEtacToapi0LhPtr(new EtacToapi0Lh(theEvtListPtr));
 MEtacToapi0Fcn fcn(theEtacToapi0LhPtr);

 MnUserParameters upar;
 if (dataSpin==2)
    {
      upar.Add("InterMass", 1.6, .1, 2.3, 0.7);
      upar.Add("InterWidth", 0.04, .01, 0.8, 0.01);
      upar.Add("spin0", 0.3, .1, 1., 0.);
      upar.Add("spin1", 0.3, .1, 1., 0.);
      upar.Add("spin2", 0.3, .1, 1., 0.);
    }
 else if (dataSpin==0)
    {
     upar.Add("InterMass", 1.1, .1, 2.3, 0.6);
     upar.Add("InterWidth", 0.04, .01, 0.8, 0.01);
     upar.Add("spin0", 0.3, .1, 1., 0.);
     upar.Add("spin1", 0.3, .1, 1., 0.);
     upar.Add("spin2", 0.3, .1, 1., 0.);
    }
 else 
   {
     Alert << "initialization of the MnUserParameters failed" << endmsg;
     exit(1);
   }

 MnMigrad migrad(fcn, upar);
 Info <<"start migrad "<< endmsg;
 FunctionMinimum min = migrad();

 if(!min.IsValid()) {
   //try with higher strategy
   Info <<"FM is invalid, try with strategy = 2."<< endmsg;
   MnMigrad migrad2(fcn, min.UserState(), MnStrategy(2));
   min = migrad2();
 }

 Info << "migrad.Fval(): " << min.Fval() << endmsg;
 
 fitParamVal theFitResult;
 theFitResult.aMass=min.UserState().Value("InterMass");
 theFitResult.aWidth=min.UserState().Value("InterWidth");
 theFitResult.cont0spin=min.UserState().Value("spin0");
 theFitResult.cont1spin=min.UserState().Value("spin1");
 theFitResult.cont2spin=min.UserState().Value("spin2");

 EtacToapi0Hist theHistogrammer(theEvtListPtr,theFitResult);

 return 0;
}

