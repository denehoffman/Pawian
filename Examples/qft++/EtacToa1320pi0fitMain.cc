#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include "Examples/qft++/EtacToa1320pi0Data.hh"
#include "Examples/qft++/EtacToa1320pi0fit.hh"
#include "Examples/qft++/EtacToa1320pi0Fcn.hh"
#include "ErrLogger/ErrLineLog.hh"
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
              << "To start the application with data containing  an intermediate resonance with Spin=0 and mass 0.98 GeV, type: ./EtacToa1320pi0fitApp -d 0\n"
              << "To start the application with data containing  an intermediate resonance with Spin=2 and mass 1.32 GeV, type: ./EtacToa1320pi0fitApp -d 2\n"
              << "with the flag -msg <errorLogMode> you can choose the mode for the error logger\n"
	      << "i.e. with './EtacToa1320pi0fitApp -d 2 -msg debugging' you start the fit with Spin=2 in the debugging mode for the error logger\n"  
	      << std::endl;
    return 0;
  }

  
  int optind=1;
  std::string dataSpinStr="";
  std::string msgModeStr="default";
  // decode arguments
  while ((optind < __argc) && (__argv[optind][0]=='-')) {

    std::string sw = __argv[optind];
    if (sw=="-d") {
      optind++;
      dataSpinStr = __argv[optind];
    }
    if (sw=="-msg"){
      optind++;
      msgModeStr = __argv[optind];
    }
    else{
      ErrMsg(warning) << "Unknown switch: " 
            << __argv[optind] << endmsg;
      optind++;
    }
  }

  ErrLineLog* myLogger=0;
  if(msgModeStr == "debugging") myLogger= new ErrLineLog(ErrLog::debugging);
  else if(msgModeStr == "trace") myLogger= new ErrLineLog(ErrLog::trace);
  else if(msgModeStr == "routine") myLogger= new ErrLineLog(ErrLog::routine);
  else if(msgModeStr == "warning")  myLogger= new ErrLineLog(ErrLog::warning);
  else if(msgModeStr == "error")    myLogger= new ErrLineLog(ErrLog::error); 
  else {
    myLogger= new ErrLineLog(ErrLog::debugging);
    ErrMsg(warning) << "ErrorLogger not (properly) set -> Use mode 'ErrLog::debugging' " << endmsg;  
  }



 std::stringstream dataSpinStrStr(dataSpinStr);
 int dataSpin=2;
 dataSpinStrStr >> dataSpin ;

 ErrMsg(routine) << "dataSpin: " << dataSpin << endmsg;

 EtacToa1320pi0fit* etacToa1320pi0fit=new EtacToa1320pi0fit(dataSpin);
 EtacToa1320pi0Fcn fcn(etacToa1320pi0fit);
 MnUserParameters upar;
 if (! etacToa1320pi0fit->initFitParameters(upar)){
   ErrMsg(fatal) << "initialization of the MnUserParameters failed in etacToa1320pi0fit->initFitParameters(upar)" << endmsg; 
 } 

 MnMigrad migrad(fcn, upar);
 ErrMsg(routine) <<"start migrad "<< endmsg;
 FunctionMinimum min = migrad();

 if(!min.IsValid()) {
   //try with higher strategy
   ErrMsg(routine) <<"FM is invalid, try with strategy = 2."<< endmsg;
   MnMigrad migrad2(fcn, min.UserState(), MnStrategy(2));
   min = migrad2();
 }

 ErrMsg(routine) << "migrad.Fval(): " << min.Fval() << endmsg;
 
//   std::cout<<"start Minos"<<std::endl;
//   MnMinos Minos(fcn, min);
//   std::pair<double,double> e0 = Minos(0);
//   std::pair<double,double> e1 = Minos(1);
//   std::pair<double,double> e2 = Minos(2);
//   std::pair<double,double> e3 = Minos(3);
//   std::pair<double,double> e4 = Minos(4);

//   std::cout<<"a1320mass: "<<min.UserState().Value("a1320mass")<<" "<<e0.first<<" "<<e0.second<<std::endl;
//   std::cout<<"a1320width: "<<min.UserState().Value("a1320width")<<" "<<e1.first<<" "<<e1.second<<std::endl;
//   std::cout<<"spin0 content: "<<min.UserState().Value("spin0")<<" "<<e2.first<<" "<<e2.second<<std::endl;
//   std::cout<<"spin1 content: "<<min.UserState().Value("spin1")<<" "<<e3.first<<" "<<e3.second<<std::endl;
//   std::cout<<"spin2 content: "<<min.UserState().Value("spin2")<<" "<<e4.first<<" "<<e4.second<<std::endl; 

  fitParamVal theFitResult;
  theFitResult.a1320Mass=min.UserState().Value("InterMass");
  theFitResult.a1320Width=min.UserState().Value("InterWidth");
  theFitResult.cont0spin=min.UserState().Value("spin0");
  theFitResult.cont1spin=min.UserState().Value("spin1");
  theFitResult.cont2spin=min.UserState().Value("spin2");
  etacToa1320pi0fit->fillFitHists(theFitResult);
  delete etacToa1320pi0fit;
  if (0!=myLogger) delete myLogger;
  return 0;
}

