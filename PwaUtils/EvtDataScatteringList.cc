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
std::string EvtDataScatteringList::DATA_PIPISCAT_NAME = "pipiScatData";
std::string EvtDataScatteringList::DATAERR_PIPISCAT_NAME = "pipiScatDataErr";
std::string EvtDataScatteringList::FIT_PIPISCAT_NAME = "pipiScatFit";
std::string EvtDataScatteringList::FITERR_PIPISCAT_NAME = "pipiScatFitErr";

EvtDataScatteringList::EvtDataScatteringList(ChannelID channelID) :
  EvtDataBaseList(channelID)
  ,_noOfInterPolData(0)
{
  if(GlobalEnv::instance()->Channel(_channelID)->channelType() != AbsChannelEnv::CHANNEL_PIPISCATTERING){
    Alert << "channelType: " << GlobalEnv::instance()->Channel(_channelID)->channelType() << " is not allowed here!!!" 
	  << "\nchannelType: " << AbsChannelEnv::CHANNEL_PIPISCATTERING << "\t needed here!!!" 
	  << endmsg;
    exit(1);
  }
 
 _noOfInterPolData=GlobalEnv::instance()->Channel(_channelID)->parser()->noOfInterpolatedScattPoints();
}


EvtDataScatteringList::~EvtDataScatteringList()
{
}

void EvtDataScatteringList::read(EventList& evtListData, EventList& evtListMc){
  if(_alreadyRead){
   Alert << "4 vectors already read " << endmsg;  // << endmsg;
    exit(1);
  }

  if(_noOfInterPolData==0){
    readScatteringData(evtListData, _evtDataList, 0);
  }
  else readScatteringDataInterpol(evtListData, _evtDataList, 0);
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
	    << "\nevtNo: " << evtData->evtNo << "\tsqrt(s): " 
	    << evtData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME)) 
	    << "\ndataPoint: " << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(DATA_PIPISCAT_NAME)) 
	    <<  "\tndataPoint error: " 
	    << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(DATAERR_PIPISCAT_NAME))
	    << "\n" << endmsg;
  }
}

void EvtDataScatteringList::readScatteringDataInterpol(EventList& evtList, 
						       std::vector<EvtData*>& theEvtList, 
						       int startNo) {
  Event* anEvent;
  int evtCount = 0;
  unsigned int noOfInterpolPoints=_noOfInterPolData;
  Event* evtOld=0;
  while ((anEvent = evtList.nextEvent())){
    if (evtCount !=0){
      std::vector<EvtData*> interpolList=convertEventInterpol(evtOld, anEvent, startNo+evtCount, noOfInterpolPoints);
      //      EvtData* evtDataInterpol=convertEventInterpol(evtOld, anEvent, startNo+evtCount);
      std::vector<EvtData*>::iterator it;
      for(it=interpolList.begin(); it!=interpolList.end(); ++it){
	theEvtList.push_back(*it);
	evtCount++;
	_noOfWeightedDataEvts++;
	InfoMsg << "added data interpolated point for pi pi scattering in PAWIAN EvtData list" 
		<< "\nevtNo: " << (*it)->evtNo << "\tsqrt(s): " 
		<< (*it)->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME)) 
		<< "\ndataPoint: " 
		<< (*it)->DoubleId.at(IdStringMapRegistry::instance()->stringId(DATA_PIPISCAT_NAME)) 
		<<  "\tndataPoint error: " 
		<< (*it)->DoubleId.at(IdStringMapRegistry::instance()->stringId(DATAERR_PIPISCAT_NAME))
		<< "\n" << endmsg;
      } 

      // theEvtList.push_back(evtDataInterpol);
      // evtCount+=noOfInterpolPoints; 
      // _noOfWeightedDataEvts++;
      // InfoMsg << "added data interpolated point for pi pi scattering in PAWIAN EvtData list" 
      // 	      << "\nevtNo: " << evtDataInterpol->evtNo << "\tsqrt(s): " << evtDataInterpol->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME)) 
      // 	      << "\ndataPoint: " << evtDataInterpol->DoubleId.at(IdStringMapRegistry::instance()->stringId(DATA_PIPISCAT_NAME)) <<  "\tndataPoint error: " << evtDataInterpol->DoubleId.at(IdStringMapRegistry::instance()->stringId(DATAERR_PIPISCAT_NAME))
      // 	      << "\n" << endmsg;
    }

    EvtData* evtData=convertEvent(anEvent, startNo+evtCount);
    theEvtList.push_back(evtData);
    ++evtCount;
    _noOfWeightedDataEvts++;
    InfoMsg << "added data point for pi pi scattering in PAWIAN EvtData list" 
	    << "\nevtNo: " << evtData->evtNo << "\tsqrt(s): " 
	    << evtData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME)) 
	    << "\ndataPoint: " 
	    << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(DATA_PIPISCAT_NAME)) 
	    <<  "\tndataPoint error: " 
	    << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(DATAERR_PIPISCAT_NAME))
	    << "\n" << endmsg;
    
    evtOld=anEvent;
  }

}

void EvtDataScatteringList::readScatteringDefaultFit(std::vector<EvtData*>& dataEvtList, 
						     std::vector<EvtData*>& fitEvtList, 
						     int startNo) {
  int evtCount = 0;
  std::vector<EvtData*>::const_iterator it=dataEvtList.begin();
  while(it!=dataEvtList.end()) {
    EvtData* evtData=new EvtData();
    evtData->evtNo=startNo+evtCount;
    evtData->DoubleMassId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME), (*it)->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME))));
    //filled now with default 0 value
    
    evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(DATA_PIPISCAT_NAME), (*it)->DoubleId.at(IdStringMapRegistry::instance()->stringId(DATA_PIPISCAT_NAME))));
    
    evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(DATAERR_PIPISCAT_NAME), (*it)->DoubleId.at(IdStringMapRegistry::instance()->stringId(DATAERR_PIPISCAT_NAME))));
    
    evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(FIT_PIPISCAT_NAME), 0.));
    
    evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(FITERR_PIPISCAT_NAME), 0.));
    
    fitEvtList.push_back(evtData); 
    
    InfoMsg << "added data point for pi pi scattering in PAWIAN EvtFit list" 
	    << "\nevtNo: " << evtData->evtNo << "\tsqrt(s): " 
	    << evtData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME)) 
	    << "\ndataPoint: " 
	    << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(DATA_PIPISCAT_NAME)) 
	    <<  "\tdataPoint error: " 
	    << evtData->DoubleId.at(IdStringMapRegistry::instance()->stringId(DATAERR_PIPISCAT_NAME))
	    << "\n" << endmsg;
    ++it;
    ++evtCount;
  }
}


EvtData* EvtDataScatteringList::convertEvent(Event* theEvent, int evtNo){
  EvtData* evtData=new EvtData();
  evtData->evtNo=evtNo;
  
  Vector4<float> current4VecFloat=*(theEvent->p4(0));
  Vector4<double> current4Vec(current4VecFloat.E(), current4VecFloat.Px(), current4VecFloat.Py(), 
			      current4VecFloat.Pz());
  evtData->DoubleMassId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME),  current4Vec.Mass()));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(DATA_PIPISCAT_NAME), theEvent->DataPoint()));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(DATAERR_PIPISCAT_NAME), theEvent->DataPointErr()));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(FIT_PIPISCAT_NAME), 0.));
  evtData->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(FITERR_PIPISCAT_NAME), 0.));
  return evtData;
}


 EvtData* EvtDataScatteringList::convertEventInterpol(Event* oldEvent, Event* newEvent, int evtNo){
 
   Vector4<float> current4VecOldFloat=*(oldEvent->p4(0));
   Vector4<float> current4VecNewFloat=*(newEvent->p4(0));
   double interpolMass=(current4VecOldFloat.Mass()+current4VecNewFloat.Mass())/2.;
   double interpolDataPoint=(oldEvent->DataPoint()+newEvent->DataPoint())/2.;
   double interpolDataErr=(oldEvent->DataPointErr()+newEvent->DataPointErr())/2.;


   EvtData* evtDataInterpol=new EvtData();   
   evtDataInterpol->evtNo=evtNo;

   evtDataInterpol->DoubleMassId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME), interpolMass));
   evtDataInterpol->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(DATA_PIPISCAT_NAME), interpolDataPoint));
   evtDataInterpol->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(DATAERR_PIPISCAT_NAME), interpolDataErr));
  evtDataInterpol->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(FIT_PIPISCAT_NAME), 0.));
evtDataInterpol->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(FITERR_PIPISCAT_NAME), 0.));


  return evtDataInterpol;
}

std::vector<EvtData*> EvtDataScatteringList::convertEventInterpol(Event* oldEvent, Event* newEvent, int evtNo, unsigned int noOfPoints){
  std::vector<EvtData*> result;
  int currentEvtNo=evtNo;

  Vector4<float> current4VecOldFloat=*(oldEvent->p4(0));
  Vector4<float> current4VecNewFloat=*(newEvent->p4(0)); 
  
  double massOld=current4VecOldFloat.Mass();
  double dataPointOld=oldEvent->DataPoint();
  double dataErrOld=oldEvent->DataPointErr();
  
  double massNew=current4VecNewFloat.Mass();
  double dataPointNew=newEvent->DataPoint();
  double dataErrNew=newEvent->DataPointErr();
  
  double slopeDataPoint=(dataPointNew-dataPointOld)/(massNew-massOld);
  double dataErrCurrent=dataErrOld+dataErrNew;
  
  double stepSize=(massNew-massOld)/(noOfPoints+1);
  
  for(double currentmass=massOld+stepSize; currentmass<massNew; currentmass+=stepSize){
    double currentDataPoint=dataPointOld-slopeDataPoint*massOld+slopeDataPoint*currentmass;
    EvtData* evtDataInterpol=new EvtData();
    evtDataInterpol->evtNo=currentEvtNo;
    evtDataInterpol->DoubleMassId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(M_PIPISCAT_NAME), currentmass));
    evtDataInterpol->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(DATA_PIPISCAT_NAME), currentDataPoint));
    evtDataInterpol->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(DATAERR_PIPISCAT_NAME), dataErrCurrent));
    evtDataInterpol->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(FIT_PIPISCAT_NAME), 0.));
    evtDataInterpol->DoubleId.insert(mapShortDouble::value_type(IdStringMapRegistry::instance()->stringId(FITERR_PIPISCAT_NAME), 0.));
    result.push_back(evtDataInterpol);
    currentEvtNo++;      
  }
  return result;  
}

