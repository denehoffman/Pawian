//************************************************************************//
//									  //
//  Copyright 2017 Markus Kuhlmann (mkuhlmann@ep1.rub.de)		  //
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
#include <sstream>
#include <vector>
#include <getopt.h>

#include <memory>

#include "PwaDynamics/BBUnstableParFit.hh"
#include "PwaDynamics/BBUnstableParFitFcn.hh"
#include "ErrLogger/ErrLogger.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

using namespace ROOT::Minuit2;

int main(int __argc,char *__argv[]){
  ErrLogger::instance().setThreshold(logging::log_level::INFO);
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){

    InfoMsg << "\nThis is a test application for calculating M and F values for BBUnstable Phasespace Calculation\n"
	    << "The switches are:\n\n"
	    << "-mRes (parameter mRes of test distribution; default 0.77)\n\n" 
	    << "-wRes (parameter wRes of test distribution; default 0.14)\n\n"
	    << "-m1 (parameter m1 of test distribution; default 0.14)\n\n" 
	    << "-m2 (parameter m2 of test distribution; default 0.14)\n\n"
	    << "-cmName (set name of the Chew-Mandelstam function; default Dudek; alternative Reid)\n\n"      
	    << endmsg;
    return 0;
  }

  std::cout << std::endl;
  std::cout << "___Starting..." << std::endl;

  // Set default values
  std::string p0Str="0.77";
  std::string p1Str="0.14";
  std::string p2Str="0.14";
  std::string p3Str="0.14";
  std::string cmName="Dudek";
  
  // Read arguments and replace default values
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-mRes"){
      optind++;
      p0Str = __argv[optind];
      found=true;
    }
    if (sw=="-wRes"){
      optind++;
      p1Str = __argv[optind];
      found=true;
    }
    if (sw=="-m1"){
      optind++;
      p2Str = __argv[optind];
      found=true;
    }
    if (sw=="-m2"){
      optind++;
      p3Str = __argv[optind];
      found=true;
    }
    if (sw=="-cmName"){
      optind++;
      cmName = __argv[optind];
      found=true;
      }
    if (!found){
      WarningMsg << "Unknown switch: " << __argv[optind] << endmsg;
      optind++;
    }
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
    }

  // Convert argument strings to doubles
  std::stringstream p0StrStr(p0Str);
  double p0=0.;
  p0StrStr >> p0;

  std::stringstream p1StrStr(p1Str);
  double p1=0.;
  p1StrStr >> p1;

  std::stringstream p2StrStr(p2Str);
  double p2=0.;
  p2StrStr >> p2;

  std::stringstream p3StrStr(p3Str);
  double p3=0.;
  p3StrStr >> p3;

  InfoMsg << "fitted with following parameter: " << endmsg;
  InfoMsg << "mRes: " << p0 << endmsg;
  InfoMsg << "wRes: " << p1 << endmsg;
  InfoMsg << "m1: " << p2 << endmsg;
  InfoMsg << "m2: " << p3 << endmsg;
  InfoMsg << "CM function name: " << cmName << endmsg;
  

  // Generate data distribution
  std::shared_ptr<BBUnstableParFit> minuitFit(new BBUnstableParFit(p0, p1, p2, p3, cmName));

  BBUnstableParFitFcn minuitFitFcn(minuitFit);

  // Set user parameters for MinuitFitFcn
  MnUserParameters upar;
  upar.Add("M", p0, 0.3*p0, 0.3, 3.0*p0);
  upar.Add("F", p1, 0.3*p1, 0.0, 5.0*p1);

  MnMigrad migrad(minuitFitFcn, upar);
  std::cout << std::endl;
  InfoMsg <<"Start Migrad "<< endmsg;
  FunctionMinimum min = migrad();
  
  if(!min.IsValid()) {
    // Try with higher strategy
    InfoMsg <<"FM is invalid, try with strategy = 2."<< endmsg;
    MnMigrad migrad2(minuitFitFcn, min.UserState(), MnStrategy(2));
    min = migrad2();
  }
  
  InfoMsg << "Function value of minimization: " << min.Fval() << endmsg;
  
  // Save final fit parameters and their errors in variables
  double final_M = min.UserState().Value("M");
  double final_F = min.UserState().Value("F");
  double eM = min.UserState().Error("M");
  double eF = min.UserState().Error("F");
  std::cout << std::endl;
  InfoMsg << "final M:\t" << final_M << " +- " << eM << endmsg;
  InfoMsg << "final F:\t" << final_F << " +- " << eF << endmsg;
  InfoMsg << "final MSqr:\t" << final_M*final_M << endmsg;
  InfoMsg << "final FSqr:\t" << final_F*final_F << endmsg;
  
  std::cout << std::endl;
  std::cout << "___Finished" << std::endl;
  std::cout << std::endl;

  return 0;

}

