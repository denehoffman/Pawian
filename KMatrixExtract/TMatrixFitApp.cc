//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
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
#include <complex>
#include <map>
#include <getopt.h>
#include <iomanip>      // std::setprecision
#include <memory>

#include "KMatrixExtract/TMatrixExtrFit.hh"
#include "KMatrixExtract/TMatrixExtrFcn.hh"

#include "ErrLogger/ErrLogger.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

using namespace ROOT::Minuit2;

int main(int __argc,char *__argv[]){
  // ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
  
  if( __argc>1 && ( strcmp( __argv[1], "-h" ) == 0 ||
                    strcmp( __argv[1], "--help" ) == 0 ) ){
    InfoMsg << "USAGE:" << endmsg;
    InfoMsg << "--path: path to kmatrix config file" << endmsg;
    InfoMsg << "--fitparams: path to fit parameter file" << endmsg;
    InfoMsg << "--sheet: sheet on which the minimization shoud be done (p corresponsds to + and n to -)" << endmsg;
    InfoMsg << "--maxImagMass: max imaginary part of the mass" << endmsg;
    InfoMsg << "--maxRealMass: max real part of the mass" << endmsg;
    InfoMsg << "--minImagMass: min imaginary part of the mass" << endmsg;
    InfoMsg << "--minRealMass: min reak part of the mass" << endmsg;
    InfoMsg << "--startRealMass: start real mass value for the minimization" << endmsg;
    InfoMsg << "--startImagMass: start imag mass value for the minimization" << endmsg;
    return 0;
  }

  std::string pathToConfigParser;
  std::string pathToFitParams="";

  std::vector<double> energyPlaneBorders;
  energyPlaneBorders.resize(4);
  energyPlaneBorders[0] = 0;    // Re min
  energyPlaneBorders[1] = 0.2; // Im min
  energyPlaneBorders[2] = 0;    // Re max
  energyPlaneBorders[3] = 0.0001;    // Im max

  std::vector<double> energyPlaneStartValue;
  energyPlaneStartValue.resize(2);
  energyPlaneStartValue[0]=0.;
  energyPlaneStartValue[1]=-0.1;  

  std::string sheetStr="";
  
  for(;optind < (__argc-1); optind++){

     std::string ws = __argv[optind];

     if(ws[0]!='-' && ws[1]!='-'){
       continue;
     }     
     else if (ws == "--path"){
	pathToConfigParser = __argv[optind+1];
     }
     else if (ws == "--fitparams"){
	pathToFitParams = __argv[optind+1];
     }
     else if (ws == "--sheet"){
       std::istringstream stream(__argv[optind+1]);
       stream >> sheetStr; 
     }
     else if(ws == "--minImagMass"){
	std::istringstream stream(__argv[optind+1]);
	stream >> energyPlaneBorders[1];
     }
     else if(ws == "--maxImagMass"){
	std::istringstream stream(__argv[optind+1]);
	stream >> energyPlaneBorders[3];
     }
     else if(ws == "--minRealMass"){
	std::istringstream stream(__argv[optind+1]);
	stream >> energyPlaneBorders[0];
     }
     else if(ws == "--maxRealMass"){
	std::istringstream stream(__argv[optind+1]);
	stream >> energyPlaneBorders[2];
     }
     else if(ws == "--startRealMass"){
	std::istringstream stream(__argv[optind+1]);
	stream >> energyPlaneStartValue[0];
     }
     else if(ws == "--startImagMass"){
	std::istringstream stream(__argv[optind+1]);
	stream >> energyPlaneStartValue[1];
     }
     else{
	Alert << "Unknown switch: " << __argv[optind] << endmsg;
	return 0;
     }
  }

  std::complex<double> energyMin(energyPlaneBorders[0], -energyPlaneBorders[1]);
  std::complex<double> energyMax(energyPlaneBorders[2], -energyPlaneBorders[3]);
  std::complex<double> energyStart(energyPlaneStartValue[0], -energyPlaneStartValue[1]);

  InfoMsg << "K-matrix config path: " << pathToConfigParser << endmsg;
  InfoMsg << "path to fitparams: " << pathToFitParams << endmsg;
  InfoMsg << "relevant sheet: " << sheetStr  << endmsg;
  InfoMsg << "energyMin: " << energyMin << endmsg;
  InfoMsg << "energyMax: " << energyMax << endmsg;
  InfoMsg << "energyStart: " << energyStart << endmsg;
  if(energyStart.real() < energyMin.real() || energyStart.real() > energyMax.real()){
    InfoMsg << "energyStart.real() is not inside the energy plane!!!" << endmsg;
    double newVal= (energyMax.real()+energyMin.real())/2.;
    energyStart= std::complex<double>(newVal, energyStart.imag()); 
    InfoMsg << "switched to energyStart: " << energyStart << endmsg;    
  }

  if(energyStart.imag() < energyMin.imag() || energyStart.imag() > energyMax.imag()){
    InfoMsg << "energyStart.imag() is not inside the energy plane!!!" << endmsg;
    double newVal= (energyMax.imag()+energyMin.imag())/2.;
    energyStart= std::complex<double>(energyStart.real(), newVal); 
    InfoMsg << "switched to energyStart: " << energyStart << endmsg;
  }

  std::shared_ptr<TMatrixExtrFit> tMatFit(new TMatrixExtrFit(pathToConfigParser, pathToFitParams, energyMin, energyMax, energyStart));
  TMatrixExtrFcn fitFcn(tMatFit);
  
  // Set user parameters for MinuitFitFcn
  MnUserParameters upar;
  upar.Add("eReal", energyStart.real(), 0.001, energyMin.real(), energyMax.real());
  upar.Add("eImag", energyStart.imag(), 0.001, energyMin.imag(), energyMax.imag());

  MnMigrad migrad(fitFcn, upar);
  InfoMsg <<"Start Migrad "<< endmsg;
  FunctionMinimum min = migrad();

  // MnMigrad migrad2(fitFcn, min.UserState(), MnStrategy(2));
  // min = migrad2();

  // MnMigrad migrad3(fitFcn, min.UserState(), MnStrategy(3));
  // min = migrad3();

  if(!min.IsValid()) {
    // Try with higher strategy
    InfoMsg <<"FM is invalid, try with strategy = 2."<< endmsg;
    MnMigrad migrad2(fitFcn, min.UserState(), MnStrategy(2));
    min = migrad2();
  }

  // Save final fit parameters and their errors in variables
  double final_eReal = min.UserState().Value("eReal");
  double final_eImag = min.UserState().Value("eImag");

  double errEReal = min.UserState().Error("eReal");
  double errEImag = min.UserState().Error("eImag");

  InfoMsg << "\n\n********************** Final fit parameters *************************\n";
  InfoMsg << std::setprecision(16) << "final eReal:\t" << final_eReal << " +- " << errEReal << endmsg;
  InfoMsg << std::setprecision(16) << "final eImag:\t" << final_eImag << " +- " << errEImag << endmsg;

  InfoMsg << "\n\n**************** Minuit FunctionMinimum information ******************" << endmsg;
    if(min.IsValid())             InfoMsg << "\n Function minimum is valid.\n";
    else                          InfoMsg << "\n WARNING: Function minimum is invalid!" << endmsg;
    if(min.HasValidCovariance())  InfoMsg << "\n Covariance matrix is valid." << endmsg;
    else                          InfoMsg << "\n WARNING: Covariance matrix is invalid!" << endmsg;
    InfoMsg <<" # of function calls: " << min.NFcn() << endmsg;
    InfoMsg <<" minimum edm: " << std::setprecision(10) << min.Edm()<<endmsg;
    if(!min.HasValidParameters()) InfoMsg << " hasValidParameters() returned FALSE" << endmsg;
    if(!min.HasAccurateCovar())   InfoMsg << " hasAccurateCovar() returned FALSE" << endmsg;
    if(!min.HasPosDefCovar()){    InfoMsg << " hasPosDefCovar() returned FALSE" << endmsg;
                                  if(min.HasMadePosDefCovar()) InfoMsg << " hasMadePosDefCovar() returned TRUE" << endmsg;
    }
    if(!min.HasCovariance())      InfoMsg << " hasCovariance() returned FALSE" << endmsg;
    if(min.HasReachedCallLimit()) InfoMsg << " hasReachedCallLimit() returned TRUE" << endmsg;
    if(min.IsAboveMaxEdm())       InfoMsg << " isAboveMaxEdm() returned TRUE" << endmsg;
    if(min.HesseFailed())         InfoMsg << " hesseFailed() returned TRUE\n" << endmsg;
  return 0;
}

