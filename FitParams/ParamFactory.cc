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

// ParamFactory class definition file. -*- C++ -*-
// Copyright 2015 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "FitParams/ParamFactory.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "FitParams/PawianParameters.hh"
#include "FitParams/MnPawianParameters.hh"

#include "ErrLogger/ErrLogger.hh"


ParamFactory* ParamFactory::_instance=0;

ParamFactory* ParamFactory::instance()
{
  if (0==_instance) _instance = new ParamFactory();
  return _instance;
}

ParamFactory::ParamFactory()
{
}

ParamFactory::~ParamFactory()
{
}

std::shared_ptr<AbsPawianParameters> ParamFactory::getParametersPointer(std::string type){
  std::shared_ptr<AbsPawianParameters> result;
  if(type=="Minuit2") result= std::shared_ptr<AbsPawianParameters>(new MnPawianParameters());
  else if(type=="Pawian") result= std::shared_ptr<AbsPawianParameters>(new PawianParameters());
  else{
    Alert << "type:\t" << type << "\tdoes not exist" << endmsg;
    exit(1);
  }
  return result;
}


