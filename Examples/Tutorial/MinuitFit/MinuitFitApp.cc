//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
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
#include <map>
#include <getopt.h>

#include <memory>

#include "Examples/Tutorial/MinuitFit/MinuitFit.hh"
#include "Examples/Tutorial/MinuitFit/MinuitFitFcn.hh"
#include "ErrLogger/ErrLogger.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

using namespace ROOT::Minuit2;

int main(int __argc,char *__argv[]){
  ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
                    || strcmp( __argv[1], "--help" ) == 0 ) ){

    InfoMsg << "\nThis is a test application for histogramming the Argand plot, the phase shift and the mass shape of the relativistic Breit-Wigner function with Blatt-Weisskopf barrier factors\n"
	 << "The switches are:\n\n"
         << "-p0 (parameter p0 of test distribution; default -10.)\n\n" 
         << "-p1 (parameter p1 of test distribution; default 10.)\n\n"
         << "-p2 (parameter p2 of test distribution; default 1.)\n\n" 
         << "-p3 (parameter p3 of test distribution; default -0.01)\n\n"
         << "-sigma (sigma of gaussian smearing; default 3)\n\n"
         << endmsg;
    return 0;
  }

  cout << endl;
  cout << "___Starting..." << endl;

  // Set default values
  std::string p0Str="-10.";
  std::string p1Str="10.";
  std::string p2Str="1.";
  std::string p3Str="-0.01";
  std::string sigmaStr="3";
  
  // Read arguments and replace default values
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-p0"){
      optind++;
      p0Str = __argv[optind];
      found=true;
    }
    if (sw=="-p1"){
      optind++;
      p1Str = __argv[optind];
      found=true;
    }
    if (sw=="-p2"){
      optind++;
      p2Str = __argv[optind];
      found=true;
    }
    if (sw=="-p3"){
      optind++;
      p3Str = __argv[optind];
      found=true;
    }
    if (sw=="-sigma"){
      optind++;
      sigmaStr = __argv[optind];
      found=true;
    }

    if (!found){
      WarningMsg << "Unknown switch: " 
            << __argv[optind] << endmsg;
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

  std::stringstream sigmaStrStr(sigmaStr);
  double sigma=0.;
  sigmaStrStr >> sigma;

  // Generate data distribution
  std::shared_ptr<MinuitFit> minuitFit(new MinuitFit(p0, p1, p2, p3, sigma));

  MinuitFitFcn minuitFitFcn(minuitFit);

  // Set user parameters for MinuitFitFcn
  MnUserParameters upar;
  upar.Add("a", -11, 3, 0., -20.);
  upar.Add("b", 9.8, 2., 15., 5.);
  upar.Add("c", 1.1, 0.3, 1.5, 0.5);
  upar.Add("d", -0.008, 0.005, 0., -0.02);

  MnMigrad migrad(minuitFitFcn, upar);
  cout << endl;
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
  double final_a = min.UserState().Value("a");
  double final_b = min.UserState().Value("b");
  double final_c = min.UserState().Value("c");
  double final_d = min.UserState().Value("d");
  double ea = min.UserState().Error("a");
  double eb = min.UserState().Error("b");
  double ec = min.UserState().Error("c");
  double ed = min.UserState().Error("d");
  cout << endl;
  InfoMsg << "final a:\t" << final_a << " +- " << ea << endmsg;
  InfoMsg << "final b:\t" << final_b << " +- " << eb << endmsg; 
  InfoMsg << "final c:\t" << final_c << " +- " << ec << endmsg;
  InfoMsg << "final d:\t" << final_d << " +- " << ed << endmsg; 
  
  // Draw the graph and the fit in one canvas
  minuitFit->drawGraph(final_a, final_b, final_c, final_d);

  cout << endl;
  cout << "___Finished" << endl;
  cout << endl;

  return 0;

}

