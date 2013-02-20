// EvtDataBaseList class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>


#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsEnv.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"

#include "Event/EventList.hh"
#include "Event/Event.hh"

#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"

#include "ErrLogger/ErrLogger.hh"


EvtDataBaseList::EvtDataBaseList(AbsEnv* theEnv) :
  _noOfWeightedDataEvts(0.),
  _noOfWeightedMcEvts(0.),
  _mcToDataRatio(1000),
  _alreadyRead(false),
  _evtNoAll(0),
  _absEnv(theEnv)
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
  read4Vecs(evtListData, _evtDataList, _noOfWeightedDataEvts, evtListData.size());

  int maxMcEvts=evtListMc.size();
  if (maxMcEvts > _mcToDataRatio*evtListData.size() ) maxMcEvts=_mcToDataRatio*evtListData.size();
  read4Vecs(evtListMc, _mcDataList, _noOfWeightedMcEvts, maxMcEvts);
  _alreadyRead=true;
}

void EvtDataBaseList::read4Vecs(EventList& evtList, std::vector<EvtData*>& theEvtList, double& evtWeightSum, int maxEvts){
  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = evtList.nextEvent())){
    if (evtCount>= maxEvts) break;
    if (evtCount%10000 == 0) Info << "4vec calculation for event " << evtCount ;  // << endmsg;

    Vector4<double> V4_all_lab(0.,0.,0.,0.);
    
    std::vector< Vector4<double> > finalState4Vecs;
    std::vector<Particle*>  finalStateParticles=_absEnv->finalStateParticles();
    std::map<std::string, Vector4<double> > particle4VecMap;
    std::vector<Particle*>::iterator itPart;
    int counter=0;
    for (itPart=finalStateParticles.begin(); itPart != finalStateParticles.end(); ++itPart){
      Vector4<float> current4VecFloat=*(anEvent->p4(counter));
      Vector4<double> current4Vec(current4VecFloat.E(), current4VecFloat.Px(), current4VecFloat.Py(), current4VecFloat.Pz());
      finalState4Vecs.push_back(current4Vec);
      particle4VecMap.insert(std::map<std::string, Vector4<double> >::value_type((*itPart)->name(), current4Vec));
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

    std::vector<Particle*>  sortedFinalStateParticles=_absEnv->finalStateParticles();    
    pawian::Collection::PtrLess thePtrLess; 
    std::sort(sortedFinalStateParticles.begin(), sortedFinalStateParticles.end(), thePtrLess);    

    std::string name_all_lab=getName(sortedFinalStateParticles);


    EvtData* evtData=new EvtData();
    evtData->FourVecsString.insert(mapString4Vec::value_type("all",V4_all_lab)); 

    //cache 4 vectors of inital state particles
    std::map<std::string, Vector4<double> >::iterator it4VecMap;
    for (it4VecMap=particle4VecMap.begin(); it4VecMap!=particle4VecMap.end(); ++it4VecMap){
      evtData->FourVecsString.insert(mapString4Vec::value_type(it4VecMap->first, it4VecMap->second));
    }
 
   //fill WignerD functions
    std::vector<boost::shared_ptr<AbsDecay> > theDecays=_absEnv->prodDecayList()->getList();
    std::vector<boost::shared_ptr<AbsDecay> >::iterator itIso;
    for (itIso=theDecays.begin(); itIso!=theDecays.end(); ++itIso){
      (*itIso)->fillWignerDs(particle4VecMap, evtData);
    }

    evtData->evtWeight=anEvent->Weight();
    evtData->evtNo=_evtNoAll;
    theEvtList.push_back(evtData);
    
    evtWeightSum += anEvent->Weight();    
    ++evtCount;
    ++_evtNoAll;
  }
}

std::string EvtDataBaseList::getName(std::vector<Particle*>& theVec){
  std::string result;

  std::vector<Particle*>::iterator it;
  for(it=theVec.begin(); it!=theVec.end(); ++it){
    result+=(*it)->name();
  }

  return result;
}
