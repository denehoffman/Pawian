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

// EvtDataBaseList class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>


#include "PwaUtils/EvtDataBaseList.hh"
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


EvtDataBaseList::EvtDataBaseList(ChannelID channelID) :
  _channelID(channelID)
  ,_noOfWeightedDataEvts(0.)
  ,_noOfWeightedMcEvts(0.)
  ,_alreadyRead(false)
  ,_finalStateParticles(GlobalEnv::instance()->Channel(_channelID)->finalStateParticles())
{
}


EvtDataBaseList::~EvtDataBaseList()
{
}

void EvtDataBaseList::read(EventList& evtListData, EventList& evtListMc){
  if(_alreadyRead){
   Alert << "4 vectors already read " << endmsg;  // << endmsg;
    exit(1);
  }
  read4Vecs(evtListData, _evtDataList, _noOfWeightedDataEvts, evtListData.size(), 0 );

  read4Vecs(evtListMc, _mcDataList, _noOfWeightedMcEvts, evtListMc.size(), evtListData.size() );
  _alreadyRead=true;
}

void EvtDataBaseList::read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList, double& evtWeightSum, int maxEvts, int startNo){
  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = evtList.nextEvent())){
    if (evtCount>= maxEvts) break;
    if (evtCount%500 == 0) Info << "4vec calculation for event " << evtCount ;  // << endmsg;

    EvtData* currentEvt=convertEvent(anEvent, startNo+evtCount);

    if (evtCount%10000 == 0){
      std::string stringAll="all";
      Vector4<double> V4_all_lab=currentEvt->FourVecsId.at(IdStringMapRegistry::instance()->stringId(stringAll));
      Info << "4vec all in lab system" << "\n"
           << " px: " << V4_all_lab.Px() <<"\t"
           << " py: " << V4_all_lab.Py() <<"\t"
           << " pz: " << V4_all_lab.Pz() <<"\t"
           << " e : " << V4_all_lab.E() << "\t"
           << " m : " << V4_all_lab.M() ;  // << endmsg;
    }

    theEvtList.push_back(currentEvt);

    evtWeightSum += anEvent->Weight();
    ++evtCount;
  }
}

EvtData* EvtDataBaseList::convertEvent(Event* theEvent, int evtNo){
    Vector4<double> V4_all_lab(0.,0.,0.,0.);

    std::vector< Vector4<double> > finalState4Vecs;
    std::map<std::string, Vector4<double> > particle4VecMap;

    std::vector<Particle*>::iterator itPart;
    int counter=0;
    for (itPart=_finalStateParticles.begin(); itPart != _finalStateParticles.end(); ++itPart){
      Vector4<float> current4VecFloat=*(theEvent->p4(counter));
      Vector4<double> current4Vec(current4VecFloat.E(), current4VecFloat.Px(), current4VecFloat.Py(), current4VecFloat.Pz());
      if(*(*itPart)==*(GlobalEnv::instance()->particleTable()->particle("photon"))){ //set mass=0
	current4Vec.SetP4(current4Vec.P(), current4Vec.Px(), current4Vec.Py(), current4Vec.Pz()); 
	}

	//      if( (current4Vec.M() != current4Vec.M()) || current4Vec.M()<1e-4) current4Vec.SetP4(sqrt(current4VecFloat.Px()*current4VecFloat.Px()+current4VecFloat.Py()*current4VecFloat.Py()+current4VecFloat.Pz()*current4VecFloat.Pz()), current4VecFloat.Px(), current4VecFloat.Py(), current4VecFloat.Pz()); 
      finalState4Vecs.push_back(current4Vec);
      particle4VecMap.insert(std::map<std::string, Vector4<double> >::value_type((*itPart)->name(), current4Vec));
      V4_all_lab += current4Vec;
      counter++;
    }

    EvtData* evtData=new EvtData();
    evtData->evtNo=evtNo;
    evtData->evtWeight=theEvent->Weight();

    std::string allString="all";
    unsigned short fourVecIdAll = IdStringMapRegistry::instance()->stringId(allString); 

    evtData->FourVecsId.insert(mapShort4Vec::value_type(fourVecIdAll,V4_all_lab));

    //cache 4 vectors of inital state particles
    std::map<std::string, Vector4<double> >::iterator it4VecMap;
    for (it4VecMap=particle4VecMap.begin(); it4VecMap!=particle4VecMap.end(); ++it4VecMap){
      std::string currentName=it4VecMap->first;
       evtData->FourVecsId.insert(mapShort4Vec::value_type(IdStringMapRegistry::instance()->stringId(currentName),  it4VecMap->second));
    }

   //fill WignerD functions
    const std::shared_ptr<AbsChannelEnv> absChannelEnvPtr=GlobalEnv::instance()->Channel(_channelID);
    std::vector<std::shared_ptr<AbsDecay> > theDecays=absChannelEnvPtr->prodDecayList()->getList();
    std::vector<std::shared_ptr<AbsDecay> >::iterator itIso;
    for (itIso=theDecays.begin(); itIso!=theDecays.end(); ++itIso){
      (*itIso)->fillWignerDs(particle4VecMap, V4_all_lab, evtData);
    }

   //fill 4Vecs for dynamics
    std::vector<std::shared_ptr<AbsDynamics> > theDynVec=DynRegistry::instance()->getDynVec();
    std::vector<std::shared_ptr<AbsDynamics> >::iterator itDyn;
    for (itDyn=theDynVec.begin(); itDyn!=theDynVec.end(); ++itDyn){
      (*itDyn)->fillMasses(evtData);
    };

   //reset filled map
    for (itIso=theDecays.begin(); itIso!=theDecays.end(); ++itIso){
      (*itIso)->resetFilledMap();
    }

    return evtData;
}

double EvtDataBaseList::noOfWeightedEvts(EventList& evtList, ChannelID channelID, int maxEvts, int startNo){
  double result=0.;
  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = evtList.nextEvent())){
    if (evtCount>= maxEvts) break;
    if (evtCount%500 == 0) Info << "4vec calculation for event " << evtCount ;  // << endmsg;

    Vector4<double> V4_all_lab(0.,0.,0.,0.);

    std::vector<Particle*>  finalStateParticles=GlobalEnv::instance()->Channel(channelID)->finalStateParticles();

    std::vector<Particle*>::iterator itPart;
    int counter=0;
    for (itPart=finalStateParticles.begin(); itPart != finalStateParticles.end(); ++itPart){
      Vector4<float> current4VecFloat=*(anEvent->p4(counter));
      Vector4<double> current4Vec(current4VecFloat.E(), current4VecFloat.Px(), current4VecFloat.Py(), current4VecFloat.Pz());
      V4_all_lab += current4Vec;
      counter++;
    }

    if (evtCount%10000 == 0){
      Info << "4vec all in lab system" << "\n"
           << " px: " << V4_all_lab.Px() <<"\t"
           << " py: " << V4_all_lab.Py() <<"\t"
           << " pz: " << V4_all_lab.Pz() <<"\t"
           << " e : " << V4_all_lab.E() << "\t"
           << " m : " << V4_all_lab.M() <<"\n"
	   << "weight: " <<  anEvent->Weight();  // << endmsg;
    }

    result += anEvent->Weight();
    ++evtCount;
  }
  return result;
}

std::string EvtDataBaseList::getName(std::vector<Particle*>& theVec){
  std::string result;

  std::vector<Particle*>::iterator it;
  for(it=theVec.begin(); it!=theVec.end(); ++it){
    result+=(*it)->name();
  }

  return result;
}
