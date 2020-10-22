//************************************************************************//
//									  //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// OmnesMatrixFactory class definition file. -*- C++ -*-
// Copyright 2016 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "qft++Extension/OmnesMatrixFactory.hh"
#include "qft++Extension/AbsOmnesMatrix.hh"
#include "qft++Extension/OmnesMatrixPiPiS5Channel.hh"
#include "qft++Extension/OmnesMatrixPiPiS3Channelwrhorho.hh"
#include "ErrLogger/ErrLogger.hh"


OmnesMatrixFactory* OmnesMatrixFactory::_instance=0;

OmnesMatrixFactory* OmnesMatrixFactory::instance()
{
  if (0==_instance) _instance = new OmnesMatrixFactory();
  return _instance;
}

OmnesMatrixFactory::OmnesMatrixFactory()
{
}

OmnesMatrixFactory::~OmnesMatrixFactory()
{
}

std::shared_ptr<AbsOmnesMatrix> OmnesMatrixFactory::getOmnesMatrixPointer(std::string type){
  std::shared_ptr<AbsOmnesMatrix> result;

  if(type=="PiPiS5Channel"){
	 result= std::shared_ptr<AbsOmnesMatrix>(new OmnesMatrixPiPiS5Channel());
 }
 else if(type=="PiPiS3Channelwrhorho"){
         result= std::shared_ptr<AbsOmnesMatrix>(new OmnesMatrixPiPiS3Channelwrhorho());
 }
 else{
    Alert << "Onmes matrix of the type:\t" << type << "\tdoes not exist" << endmsg;
    exit(1);
  }

  return result;
}

