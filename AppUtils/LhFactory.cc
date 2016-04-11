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

// LhFactory class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "AppUtils/LhFactory.hh"

#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/GlobalEnv.hh"

#include "pbarpUtils/pbarpHeliLh.hh"
#include "pbarpUtils/pbarpCanoLh.hh"
#include "pbarpUtils/pbarpTensorLh.hh"
#include "epemUtils/epemHeliLh.hh"
#include "epemUtils/epemTensorLh.hh"
#include "resUtils/resBaseLh.hh"
#include "ErrLogger/ErrLogger.hh"


LhFactory* LhFactory::_instance=0;

LhFactory* LhFactory::instance()
{
  if (0==_instance) _instance = new LhFactory();
  return _instance;
}

LhFactory::LhFactory()
{
}

LhFactory::~LhFactory()
{
}

std::shared_ptr<AbsLh> LhFactory::getLh(short channelType, ChannelID id, std::string& formalism){
  std::shared_ptr<AbsLh> result;  
  if( channelType == AbsChannelEnv::CHANNEL_PBARP){
    if(formalism =="Cano") result=std::shared_ptr<AbsLh>(new pbarpCanoLh(id));
    else if(formalism == "Heli") result=std::shared_ptr<AbsLh>(new pbarpHeliLh(id));
    else if(formalism == "Tensor") result=std::shared_ptr<AbsLh>(new pbarpTensorLh(id));
    else {
      Alert << "prodFormalism\t" << formalism << "\tfor channel type AbsChannelEnv::CHANNEL_PBARP doesn't exist!!!" << endmsg;
      exit(1);
    }
  }
  else if( channelType == AbsChannelEnv::CHANNEL_EPEM){
    if (formalism=="Heli" || formalism=="HeliMultipole") result = std::shared_ptr<AbsLh>(new epemHeliLh(id));
    else if (formalism=="Tensor") result = std::shared_ptr<AbsLh>(new epemTensorLh(id));
    else {
      Alert << "prodFormalism\t" << formalism << "\tfor channel type AbsChannelEnv::CHANNEL_EPEM doesn't exist!!!" << endmsg;
      exit(1);
    }
  }
  else if( channelType == AbsChannelEnv::CHANNEL_RES){
    result = std::shared_ptr<AbsLh>(new resBaseLh(id));
  }
  
  return result;
}


