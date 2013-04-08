//************************************************************************//
//                                                                        //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 Julian Pychy (julian@ep1.rub.de)                       //
//                 - Ruhr-Universität Bochum                              //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

#include <getopt.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <thread>

#include <boost/thread.hpp>

#include "PwaUtils/PwaGen.hh"
#include "PwaUtils/AbsEnv.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/ParserBase.hh"

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include "PspGen/EvtGenKine.hh"
#include "PspGen/EvtRandom.hh"
#include "PspGen/EvtRandomEngine.hh"
#include "PspGen/EvtSimpleRandomEngine.hh"
#include "PspGen/EvtVector4R.hh"

#include "Event/Event.hh"
#include "Event/EventList.hh"

#include "TFile.h"
#include "TH1F.h"

PwaGen::PwaGen(AbsEnv* theEnv) :
  _absEnv(theEnv)
  ,_initial4Vec(EvtVector4R(theEnv->initial4Vec().E(), theEnv->initial4Vec().Px(), theEnv->initial4Vec().Py(), theEnv->initial4Vec().Pz()))
  ,_finalStateParticles(theEnv->finalStateParticles())
  ,_stream(new std::ofstream("evtGen.dat"))
  ,_genWithModel(true)
  ,_unitScaleFactor(1.)
  ,_energyFirst(false)
  ,_useEvtWeight(theEnv->parser()->useEvtWeight())
{
  _theTFile=new TFile("EvtGen.root","recreate");
  inv01MassH1=new TH1F("inv01MassH1","inv01MassH1",500, 0., 3.);
  inv02MassH1=new TH1F("inv02MassH1","inv02MassH1",500, 0., 3.);
  inv12MassH1=new TH1F("inv12MassH1","inv12MassH1",500, 0., 3.);
  _genWithModel=theEnv->parser()->generateWithModel();

  std::string unit=theEnv->parser()->unitInFile();

  if(unit=="GEV"){
    _unitScaleFactor=1.;
  }
  else if(unit=="MEV"){
    _unitScaleFactor=1000.;
    }
  else{
    Alert << "unit " << unit << " does not exist!!!" <<endmsg;
    exit(0);
  }

  std::string order =theEnv->parser()->orderInFile();
  if(order=="E Px Py Pz"){
    _energyFirst=true;
  }
  else if(order=="Px Py Pz E"){
    _energyFirst=false;
  }
  else{
    Alert << "order " << order << " does not exist!!!" <<endmsg;
    exit(0);
  }
}

PwaGen::~PwaGen()
{
  _theTFile->Write();
  _theTFile->Close();
  _stream->close();
}

void PwaGen::generate(boost::shared_ptr<AbsLh> theLh, fitParams& theFitParams){
  // std::ofstream theStream ("output.dat");
  int counterFsp=0;
  std::vector<Particle*>::iterator it;
  for(it=_finalStateParticles.begin(); it!=_finalStateParticles.end(); ++it){
    mass[counterFsp]=(*it)->mass();
    counterFsp++;
  }

  EvtSimpleRandomEngine myRandom(_absEnv->parser()->randomSeed());
  EvtRandom::setRandomEngine(&myRandom);
  bool generateEvents=true;
  int noOfAcceptedEvts=0;
  int noOfEvtsToGenerate=_absEnv->parser()->noOfGenEvts();
  int noOfAllGenEvts=0;
  int noOfIterations=0;

  while(generateEvents){
    noOfIterations++;

    EventList currentEvtList;
    for (int count = 0; count < 100000; count++){
      noOfAllGenEvts++;
      EvtVector4R p4[_finalStateParticles.size()];
      EvtGenKine::PhaseSpace(_finalStateParticles.size(), mass, p4, _initial4Vec.mass());
      for(unsigned int t=0; t<_finalStateParticles.size(); ++t){
	p4[t].applyBoostTo(_initial4Vec);
      }
      //    dumpAscii(p4);
      addEvt(currentEvtList, p4, count);
      
      //    EvtData* theData
      inv01MassH1->Fill((p4[0]+p4[1]).mass());
      inv02MassH1->Fill((p4[0]+p4[2]).mass());
      inv12MassH1->Fill((p4[1]+p4[2]).mass());
    }
    
    currentEvtList.rewind();
    
    Event* anEvent;
    int evtCount = 0;
    while ((anEvent = currentEvtList.nextEvent()) != 0 && evtCount < 10) {
      Info        << "\n";
      for(unsigned int i=0; i<_finalStateParticles.size(); ++i){
	Info        << (*anEvent->p4(i)) << "\tm = " << anEvent->p4(i)->Mass() << "\n";
      }
      Info        << "\n" << endmsg;
      ;  // << endmsg;
      ++evtCount;
    }
    
    currentEvtList.rewind();
    
    boost::shared_ptr<EvtDataBaseList> eventListPtr(new EvtDataBaseList(_absEnv));
    
    std::vector<EvtData*> dataList;
    double evtWeightSum=0.;
    eventListPtr->read4Vecs(currentEvtList, dataList, evtWeightSum, 100000, noOfAllGenEvts-100000); 
    
    theLh->updateFitParams(theFitParams);
    std::vector<EvtData*>::const_iterator itEvt;
 
   if(!_genWithModel){
      itEvt=dataList.begin();
      while(itEvt!=dataList.end()){
	dumpAscii(*itEvt);
	noOfAcceptedEvts++;
	if(noOfAcceptedEvts==noOfEvtsToGenerate){
	  generateEvents=false;
	  break;
	}
	++itEvt;
      }
    }
   else{ //generate with model    
     //fit retrieve maxFitWeight
     double maxFitWeight=0.;
     itEvt=dataList.begin();
     while(itEvt!=dataList.end())
       {
	 double fitWeight= theLh->calcEvtIntensity( *itEvt, theFitParams );
	 DebugMsg << (*itEvt)->evtNo <<  "\tfitWeight:\t" << fitWeight << endmsg;
	 if (maxFitWeight< fitWeight) maxFitWeight=fitWeight;
	 ++itEvt;
       }
     
     itEvt=dataList.begin();
     while(itEvt!=dataList.end())
       {
	 double fitWeight= theLh->calcEvtIntensity( *itEvt, theFitParams );
	 double randWeight = EvtRandom::Flat( 0., maxFitWeight );
	 if ( randWeight < fitWeight ){
	   dumpAscii(*itEvt);
	   noOfAcceptedEvts++;
	   if(noOfAcceptedEvts==noOfEvtsToGenerate){
	     generateEvents=false;
	     break;
	   }
	 }
	 ++itEvt;
       }
     
   }
   Info << "iteration:\t" << noOfIterations << "\tnoOfAcceptedEvts:\t" << noOfAcceptedEvts << endmsg;

   //delete content of data list
   for(itEvt=dataList.begin(); itEvt!=dataList.end(); ++itEvt) delete (*itEvt); 
  }  
  
}

void PwaGen::addEvt(EventList& evtList, EvtVector4R* evt4Vec4Rs,int evtNumber){
  Event* newEvent = new Event();
  newEvent->addWeight(1.);
  for(unsigned int t=0; t<_finalStateParticles.size(); ++t){
    newEvent->addParticle(evt4Vec4Rs[t].get(0), evt4Vec4Rs[t].get(1), evt4Vec4Rs[t].get(2), evt4Vec4Rs[t].get(3));
  }

  evtList.add(newEvent);
}


void PwaGen::dumpAscii(EvtData* evtData){
  if (_useEvtWeight) *_stream << 1.000 << std::endl;
  std::vector<Particle* > fsParticles = _absEnv->finalStateParticles();

  std::vector<Particle* >::const_iterator fspIter = fsParticles.begin();
  for( ; fspIter != fsParticles.end(); ++fspIter ) {
    Vector4<double> tmp4vec = evtData->FourVecsString[ (*fspIter)->name() ];
    if(_energyFirst){ 
      *_stream << std::setprecision(8)  << tmp4vec.E()*_unitScaleFactor  << "\t" << tmp4vec.Px()*_unitScaleFactor << "\t" << tmp4vec.Py()*_unitScaleFactor << "\t" << tmp4vec.Pz()*_unitScaleFactor << "\t" << std::endl;
    }
    else{
      *_stream << std::setprecision(8)  << tmp4vec.Px()*_unitScaleFactor << "\t" << tmp4vec.Py()*_unitScaleFactor << "\t" << tmp4vec.Pz()*_unitScaleFactor << "\t" << tmp4vec.E()*_unitScaleFactor << std::endl;
    }
  }
}
