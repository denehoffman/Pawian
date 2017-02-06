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

// EvtDataScatteringList class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>


#include "PwaUtils/EvtDataScatteringList.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/AbsDynamics.hh"
#include "PwaUtils/DynRegistry.hh"
#include "PwaUtils/GlobalEnv.hh"

#include "Event/EventList.hh"
#include "Event/Event.hh"

#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/IdStringMapRegistry.hh"

#include "ConfigParser/ParserBase.hh"

#include "ErrLogger/ErrLogger.hh"

std::string EvtDataScatteringList::M_PIPISCAT_NAME = "pipiScatm";
std::string EvtDataScatteringList::PHI_PIPISCAT_NAME= "pipiScatphi";
std::string EvtDataScatteringList::PHIERR_PIPISCAT_NAME= "pipiScatphiErr";
std::string EvtDataScatteringList::ETA_PIPISCAT_NAME= "pipiScateta";
std::string EvtDataScatteringList::ETAERR_PIPISCAT_NAME= "pipiScateatErr";
std::string EvtDataScatteringList::PHIFIT_PIPISCAT_NAME= "pipiScatphiFit";
std::string EvtDataScatteringList::ETAFIT_PIPISCAT_NAME= "pipiScatetaFit";

EvtDataScatteringList::EvtDataScatteringList(ChannelID channelID) :
  EvtDataBaseList(channelID)
{
  if(GlobalEnv::instance()->Channel(_channelID)->channelType() != AbsChannelEnv::CHANNEL_PIPISCATTERING){
    Alert << "channelType: " << GlobalEnv::instance()->Channel(_channelID)->channelType() << " is not allowed here!!!" 
	  << "\nchannelType: " << AbsChannelEnv::CHANNEL_PIPISCATTERING << "\t needed here!!!" 
	  << endmsg;
    exit(1);
  }
}


EvtDataScatteringList::~EvtDataScatteringList()
{
}

void EvtDataScatteringList::read(EventList& evtListData, EventList& evtListMc){
  if(_alreadyRead){
   Alert << "4 vectors already read " << endmsg;  // << endmsg;
    exit(1);
  }

  readScatteringData(evtListData, _evtDataList, 0);
  readScatteringDefaultFit(_evtDataList, _mcDataList, evtListData.size());
  _alreadyRead=true;
}

void EvtDataScatteringList::readScatteringData(EventList& evtList, std::vector<EvtData*>& theEvtList, int startNo){
  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = evtList.nextEvent())){
    EvtData* evtData=convertEvent(anEvent, startNo+evtCount);
    theEvtList.push_back(evtData);
    ++evtCount;
    _noOfWeightedDataEvts++;
    InfoMsg << "added data point for pi pi scattering in PAWIAN EvtData list" 
	    << "\nevtNo: " << evtData->evtNo << "\tsqrt(s): " << evtData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME)) 
	    << "\nphi: " << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(PHI_PIPISCAT_NAME)) <<  "\tphi error: " << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(PHIERR_PIPISCAT_NAME))
	    << "\neta: " << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(ETA_PIPISCAT_NAME)) <<  "\teta error: " << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(ETAERR_PIPISCAT_NAME))
	    << "\n" << endmsg;
  }
}

void EvtDataScatteringList::readScatteringDefaultFit(std::vector<EvtData*>& dataEvtList, std::vector<EvtData*>& fitEvtList, int startNo){
  int evtCount = 0;
  std::vector<EvtData*>::const_iterator it=dataEvtList.begin();
  while(it!=dataEvtList.end())
    {
      EvtData* evtData=new EvtData();
      evtData->evtNo=startNo+evtCount;
      evtData->DoubleMassId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME), (*it)->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME))));
      //filled now with default 0 value

      evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(PHI_PIPISCAT_NAME), (*it)->DoubleId.at(IdStringMapRegistry::instance()->stringId(PHI_PIPISCAT_NAME))));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(PHIERR_PIPISCAT_NAME), 0.));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(ETA_PIPISCAT_NAME), 0.));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(ETAERR_PIPISCAT_NAME), 0.));
      fitEvtList.push_back(evtData)
;    InfoMsg << "added data point for pi pi scattering in PAWIAN EvtFit list" 
	    << "\nevtNo: " << evtData->evtNo << "\tsqrt(s): " << evtData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME)) 
	    << "\nphi: " << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(PHI_PIPISCAT_NAME)) <<  "\tphi error: " << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(PHIERR_PIPISCAT_NAME))
	    << "\neta: " << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(ETA_PIPISCAT_NAME)) <<  "\teta error: " << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(ETAERR_PIPISCAT_NAME))
	    << "\n" << endmsg;
      ++it;
      ++evtCount;
    }
}


EvtData* EvtDataScatteringList::convertEvent(Event* theEvent, int evtNo){
  EvtData* evtData=new EvtData();
  evtData->evtNo=evtNo;
  
  Vector4<float> current4VecFloat=*(theEvent->p4(0));
  Vector4<double> current4Vec(current4VecFloat.E(), current4VecFloat.Px(), current4VecFloat.Py(), current4VecFloat.Pz());
  evtData->DoubleMassId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME),  current4Vec.Mass()));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(PHI_PIPISCAT_NAME), theEvent->Phase()));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(PHIERR_PIPISCAT_NAME), theEvent->PhaseErr()));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(ETA_PIPISCAT_NAME), theEvent->Eta()));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(ETAERR_PIPISCAT_NAME), theEvent->EtaErr()));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(PHIFIT_PIPISCAT_NAME), 0.));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(ETAFIT_PIPISCAT_NAME), 0.));
  return evtData;
}

