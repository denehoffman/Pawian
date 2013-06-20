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

// AppBase class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/AppBase.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/AbsEnv.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/PwaGen.hh"

#include "ErrLogger/ErrLogger.hh"


AppBase::AppBase(AbsEnv* absEnv, std::shared_ptr<AbsLh> theLhPtr, std::shared_ptr<FitParamsBase> theFitParamBase) :
  _absEnv(absEnv), 
  _absLhPtr(theLhPtr),
  _fitParamBasePtr(theFitParamBase)
{
}

AppBase::~AppBase()
{
}

void AppBase::dumpDefaultParams(){
    fitParams defaultVal;
    fitParams defaultErr;
    _absLhPtr->getDefaultParams(defaultVal, defaultErr);

    std::stringstream defaultparamsname;
    defaultparamsname << "defaultparams" << _absEnv->outputFileNameSuffix() << ".dat";
    std::ofstream theStreamDefault ( defaultparamsname.str().c_str() );
    
    _fitParamBasePtr->dumpParams(theStreamDefault, defaultVal, defaultErr);
}

void AppBase::generate(fitParams& theParams){
    std::shared_ptr<PwaGen> pwaGenPtr(new PwaGen(_absEnv));
    pwaGenPtr->generate(_absLhPtr, theParams);
    _fitParamBasePtr->printParams(theParams);
}

