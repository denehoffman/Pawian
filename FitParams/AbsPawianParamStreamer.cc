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

#include <fstream>

#include "FitParams/AbsPawianParamStreamer.hh"

#include "ErrLogger/ErrLogger.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "FitParams/ParamFactory.hh"

AbsPawianParamStreamer::AbsPawianParamStreamer(std::string& filePath) :
  _pawianParams(ParamFactory::instance()->getParametersPointer("Pawian"))
{
  //  _pawianParams=ParamFactory::instance()->getParametersPointer("Pawian");
  std::ifstream ifs(filePath.c_str());
  if(!ifs.good()) 
    {
      Alert << "Error: could not parse " << filePath << endmsg;
      exit(1);
    }

  for (std::string lineStr; std::getline(ifs, lineStr); ){
    std::string paramName, strTmp;
    double val, err, min, max;
    min=0.;
    max=-1.;
    err=0.;

    Info << lineStr << endmsg;
    istringstream lineIss(lineStr);
    
    lineIss >> paramName;
    lineIss >> strTmp;
    val=atof(strTmp.c_str());
    
    if (!lineIss.eof()){
    lineIss >> strTmp;
    err=atof(strTmp.c_str());
    }

    if (!lineIss.eof()){
      lineIss >> strTmp;
      min=atof(strTmp.c_str());
    }

    if (!lineIss.eof()){
      lineIss >> strTmp;
      max=atof(strTmp.c_str());
    }
 
    Info << "paramName: "<< paramName << "\tval: " << val << "\terr:" << err << "\tmin:" << min << "\tmax:" << max << endmsg;

    fillDefault(paramName, val, err, min, max); 
  }
  std::cout << "\n\n parameters in pawianParams are: \n\n";
  _pawianParams->print(std::cout);
}

AbsPawianParamStreamer::~AbsPawianParamStreamer(){;}

void AbsPawianParamStreamer::fillDefault(const std::string &theName, double val, double err, double min, double max){
  _pawianParams->Add(theName, val, err);
  if(min<max)  _pawianParams->SetLimits(theName, min, max); 
  return;
}

