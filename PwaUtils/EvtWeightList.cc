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

// EvWeightList class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#include <getopt.h>


#include "PwaUtils/EvtWeightList.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"

EvtWeightList::EvtWeightList(ChannelID channelID) :
  EvtDataBaseList(channelID)
{
}


EvtWeightList::~EvtWeightList()
{
}

void EvtWeightList::read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList, double& evtWeightSum, int maxEvts, int startNo){
  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = evtList.nextEvent())){
    if (evtCount>= maxEvts) break;
    //    if (evtCount%10000 == 0) Info << "4vec calculation for event " << evtCount ;  // << endmsg;
    if (evtCount%500 == 0) Info << "4vec calculation for event " << evtCount ;  // << endmsg;

    Vector4<double> V4_all_lab(0.,0.,0.,0.);

    std::vector<Particle*>  finalStateParticles=GlobalEnv::instance()->Channel(_channelID)->finalStateParticles();

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
           << " m : " << V4_all_lab.M() ;  // << endmsg;
    }

    // EvtData* evtData=new EvtData();
    // evtData->evtWeight=anEvent->Weight();
    // evtData->evtNo=startNo+evtCount;
    //    theEvtList.push_back(evtData);

    evtWeightSum += anEvent->Weight();
    ++evtCount;
  }
}

