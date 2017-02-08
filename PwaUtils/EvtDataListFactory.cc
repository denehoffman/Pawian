//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// EvtDataListFactory class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/EvtDataListFactory.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "EvtDataScatteringList.hh"
#include "ErrLogger/ErrLogger.hh"


EvtDataListFactory* EvtDataListFactory::_instance=0;

EvtDataListFactory* EvtDataListFactory::instance()
{
  if (0==_instance) _instance = new EvtDataListFactory();
  return _instance;
}

EvtDataListFactory::EvtDataListFactory()
{
}

EvtDataListFactory::~EvtDataListFactory()
{
}

std::shared_ptr<EvtDataBaseList> EvtDataListFactory::evtDataListPtr(std::shared_ptr<AbsChannelEnv> absChannelEnv){
  std::shared_ptr<EvtDataBaseList> eventListPtr;
  if(absChannelEnv->channelType() == AbsChannelEnv::CHANNEL_PIPISCATTERING) eventListPtr=std::shared_ptr<EvtDataBaseList>(new EvtDataScatteringList(absChannelEnv->channelID()));
  else eventListPtr=std::shared_ptr<EvtDataBaseList>(new EvtDataBaseList(absChannelEnv->channelID()));
  return eventListPtr; 
}
