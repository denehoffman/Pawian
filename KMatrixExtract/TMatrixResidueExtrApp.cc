//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
//                 Xiaoshuai Qin (xqin@ep1.rub.de)                        //
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
#include "KMatrixExtract/TMatrixResidueExtr.hh"
#include "KMatrixExtract/TMatrixResiduePathExtr.hh"
#include "FVectorResidueExtr.hh"
#include "KMatrixExtract/FVectorResiduePathExtr.hh"
#include "ConfigParser/pipiScatteringParser.hh"

#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
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

 std::shared_ptr<TMatrixResidueExtr> tMatrixResidueExtr;
 if (productionFormalism =="PVecIntensity"){
     //tMatrixResidueExtr = std::shared_ptr<FVectorResidueExtr>(new FVectorResidueExtr(theParser));

     if(theParser->residueExtrMethod()=="Laurent") tMatrixResidueExtr = std::shared_ptr<FVectorResidueExtr>(new FVectorResidueExtr(theParser));
     else if(theParser->residueExtrMethod()=="Cauchy") tMatrixResidueExtr = std::shared_ptr<FVectorResidueExtr>(new FVectorResiduePathExtr(theParser));
     else{
         Alert << "residue extraction method: " << theParser->residueExtrMethod() << " is not supported!!!!" << endmsg;
         exit(0);
     }

 }
  else{
    if(theParser->residueExtrMethod()=="Laurent") tMatrixResidueExtr = std::shared_ptr<TMatrixResidueExtr> (new TMatrixResidueExtr(theParser));
    else if(theParser->residueExtrMethod()=="Cauchy") tMatrixResidueExtr = std::shared_ptr<TMatrixResidueExtr> (new TMatrixResiduePathExtr(theParser));
    else{
      Alert << "residue extraction method: " << theParser->residueExtrMethod() << " is not supported!!!!" << endmsg;
      exit(0);
    }
  }

  tMatrixResidueExtr->Calculation();
  return 0;
}
