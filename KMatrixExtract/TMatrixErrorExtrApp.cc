//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
//                 Xiaoshuai Qin (xqin@ep1.rub.de)             //
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

#include "FitParams/PwaCovMatrix.hh"
#include "KMatrixExtract/TMatrixErrorExtr.hh"
#include "KMatrixExtract/TMatrixExtrFcn.hh"

#include "ErrLogger/ErrLogger.hh"
#include "ConfigParser/pipiScatteringParser.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

using namespace ROOT::Minuit2;

int main(int __argc,char *__argv[]){
  // ErrLogger::instance().setThreshold(logging::log_level::DEBUG);
  ErrLogger::instance().setThreshold(logging::log_level::INFO);


  if( __argc>1 && ( strcmp( __argv[1], "-h" ) == 0 ||
		    strcmp( __argv[1], "--help" ) == 0 ) ){
    InfoMsg << "USAGE:" << endmsg;
    InfoMsg << "-c: path to config file" << endmsg;
    InfoMsg << "additional options are also possible" << endmsg;
    return 0;
  }

  pipiScatteringParser* theParser = new pipiScatteringParser(__argc, __argv);
  InfoMsg << "pathToKMatrixCompareFile: " << theParser->pathToKMatrixCompareFile() << endmsg;

  std::string productionFormalism=theParser->productionFormalism();
  const std::vector<std::string> decayDynamicsVec=theParser->decayDynamics();
  InfoMsg << "decayDynamicsVec.size(): " << decayDynamicsVec.size() << endmsg;

  if(decayDynamicsVec.size() != 1){
    Alert << "the number of decay dynamics is " <<  decayDynamicsVec.size() << endmsg;
    Alert << "the number is required to be exatly 1!!! " << endmsg;
    return 0;
  }

  if (productionFormalism!="PVecIntensity" && productionFormalism!="Phase" && productionFormalism!="Elasticity" 
      && productionFormalism!="ArgandUnits"){
    Alert << "production formalism " << productionFormalism << " is not allowed!!!" << endmsg;
    InfoMsg << "It is required to use PVecIntensity, Phase, Elasticity or ArgandUnits !!!!" << endmsg;
    return 0;    
  }  

  std::shared_ptr<TMatrixErrorExtr> tMatErrExtract;
  if (productionFormalism =="PVecIntensity"){
    return 0;
  }
  else{
    tMatErrExtract = std::shared_ptr<TMatrixErrorExtr> (new TMatrixErrorExtr(theParser));
  }

  tMatErrExtract->Calculation();
  tMatErrExtract->printErrors();
  std::complex<double> error = tMatErrExtract->GetError();
  std::complex<double> result = tMatErrExtract->GetResult();
  InfoMsg << "      " << endmsg;
  InfoMsg << "      " << endmsg;
  InfoMsg << "      " << endmsg;
  InfoMsg << "result (M  -i/2 Gamma): " << result.real() << " +/- " << error.real() << "\t" << result.imag() << " +/- " << error.imag() << endmsg;
  InfoMsg << "result (M  Gamma): " << result.real() << " +/- " << error.real() << "\t" << -2.*result.imag() << " +/- " << 2.*error.imag() << endmsg;
  return 0;
}
