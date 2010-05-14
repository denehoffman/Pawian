#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include "Examples/qft++/EtacToa1320pi0Data.hh"
#include "Examples/qft++/EtacToa1320pi0fit.hh"
#include "Examples/qft++/EtacToa1320pi0Fcn.hh"
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
	      << std::endl;
    return 0;
  }

  int optind=1;
  std::string dataSpinStr="";
  // decode arguments
  while ((optind < __argc) && (__argv[optind][0]=='-')) {

    std::string sw = __argv[optind];
    if (sw=="-d") {
      optind++;
      dataSpinStr = __argv[optind];
    }
    else{
      cout << "Unknown switch: " 
           << __argv[optind] << endl;
      optind++;
    }
  }

 std::stringstream dataSpinStrStr(dataSpinStr);
 int dataSpin=2;
 dataSpinStrStr >> dataSpin ;

  std::cout << "dataSpin: " << dataSpin << std::endl;

 EtacToa1320pi0fit* etacToa1320pi0fit=new EtacToa1320pi0fit(dataSpin);
 EtacToa1320pi0Fcn fcn(etacToa1320pi0fit);
 MnUserParameters upar;
 if (! etacToa1320pi0fit->initFitParameters(upar)) assert(0);

 MnMigrad migrad(fcn, upar);
 std::cout<<"start migrad "<<std::endl;
 FunctionMinimum min = migrad();

 if(!min.IsValid()) {
   //try with higher strategy
   std::cout<<"FM is invalid, try with strategy = 2."<<std::endl;
   MnMigrad migrad2(fcn, min.UserState(), MnStrategy(2));
   min = migrad2();
 }

//  std::cout << "minimum: " << min << std::endl;
 std::cout << "migrad.Fval(): " << min.Fval() << std::endl;
 
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

  return 0;
}

