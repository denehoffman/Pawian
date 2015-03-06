//************************************************************************//
//                                                                        //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 Julian Pychy (julian@ep1.rub.de)                       //
//                 - Ruhr-Universit??t Bochum                              //
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

#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>

#include "PwaUtils/PwaGen.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "ConfigParser/ParserBase.hh"
#include "PwaUtils/GlobalEnv.hh"

#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include "PspGen/EvtGenKine.hh"
#include "PspGen/EvtRandom.hh"
#include "PspGen/EvtMTRandomEngine.hh"
#include "PspGen/EvtVector4R.hh"

#include "Event/Event.hh"
#include "Event/EventList.hh"
#include "Event/MassRangeCut.hh"

#include "TFile.h"
#include "TH1F.h"

const unsigned int EVENTS_PER_ITERATION = 100000;

PwaGen::PwaGen() :
   _energyFirst(false)
  ,_useEvtWeight(GlobalEnv::instance()->parser()->useEvtWeight())
  ,_genWithModel(GlobalEnv::instance()->parser()->generateWithModel())
  ,_unitScaleFactor(1.)
  , _maxFitWeight(0.)
  ,_initial4Vec(EvtVector4R(GlobalEnv::instance()->Channel()->initial4Vec().E(), GlobalEnv::instance()->Channel()->initial4Vec().Px(),
                            GlobalEnv::instance()->Channel()->initial4Vec().Py(), GlobalEnv::instance()->Channel()->initial4Vec().Pz()))
  ,_finalStateParticles(GlobalEnv::instance()->Channel()->finalStateParticles())
 {
  std::ostringstream datFileName;
  datFileName << "evtGen" << GlobalEnv::instance()->outputFileNameSuffix() << ".dat";
  _stream = new std::ofstream(datFileName.str());
  std::ostringstream rootFileName;
  rootFileName << "evtGen" << GlobalEnv::instance()->outputFileNameSuffix() << ".root";
  _theTFile=new TFile(rootFileName.str().c_str(),"recreate");

  inv01MassH1=new TH1F("inv01MassH1","inv01MassH1",500, 0., 3.);
  inv02MassH1=new TH1F("inv02MassH1","inv02MassH1",500, 0., 3.);
  inv12MassH1=new TH1F("inv12MassH1","inv12MassH1",500, 0., 3.);
  inv01MassWeightH1=new TH1F("inv01MassWeightH1","inv01MassWeightH1",500, 0., 3.);
  inv02MassWeightH1=new TH1F("inv02MassWeightH1","inv02MassWeightH1",500, 0., 3.);
  inv12MassWeightH1=new TH1F("inv12MassWeightH1","inv12MassWeightH1",500, 0., 3.);

  std::string unit=GlobalEnv::instance()->parser()->unitInFile();

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

  std::string order = GlobalEnv::instance()->parser()->orderInFile();
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

  for(unsigned int i=0; i<_finalStateParticles.size(); i++){
     _fspMasses[i]=_finalStateParticles.at(i)->mass();
  }
}

PwaGen::~PwaGen()
{
  _theTFile->Write();
  _theTFile->Close();
  _stream->close();
}




std::shared_ptr<EventList> PwaGen::GeneratePspEventList(unsigned int numEvents){
  bool useMassRange = GlobalEnv::instance()->Channel()->useMassRange();

  std::vector< std::shared_ptr<MassRangeCut> > massRangeCuts= GlobalEnv::instance()->Channel()->massRangeCuts();
  std::vector< std::shared_ptr<MassRangeCut> >::iterator itMassRangeCut;

  std::shared_ptr<EventList> eventList(new EventList);

  for(unsigned int i=0; i<numEvents; i++){
     EvtVector4R p4[_finalStateParticles.size()];
     EvtGenKine::PhaseSpace(_finalStateParticles.size(), _fspMasses, p4, _initial4Vec.mass());
     for(unsigned int j=0; j<_finalStateParticles.size(); j++){
	p4[j].applyBoostTo(_initial4Vec);
     }

     bool acceptEvt=true;
     if(useMassRange){
       for (itMassRangeCut=massRangeCuts.begin(); itMassRangeCut!=massRangeCuts.end(); ++itMassRangeCut){
	 EvtVector4R particleSystem(0,0,0,0);
	 std::vector<unsigned int> particleIndices=(*itMassRangeCut)->particleIds();
	 for(auto it = particleIndices.begin(); it!=particleIndices.end(); ++it){
	   particleSystem = particleSystem + p4[*it];
	 }
	 double invMass = particleSystem.mass();
//         Info << "invMass: " << invMass << endmsg; 
	 if(invMass < (*itMassRangeCut)->massMin() || invMass > (*itMassRangeCut)->massMax()){
	   acceptEvt=false;
//	   Info << "event not accepted" << endmsg; 
	   break;
	 }
       }
     } 

     if(acceptEvt){
       Info << "event no " << i << " accepted" << endmsg;
       AddEventToEventList(eventList, p4, i);
       
       inv01MassH1->Fill((p4[0]+p4[1]).mass());
       if(_finalStateParticles.size()>2){
	 inv02MassH1->Fill((p4[0]+p4[2]).mass());
	 inv12MassH1->Fill((p4[1]+p4[2]).mass());
       }
     }
     else{
       i--;
     }
  }
     
    eventList->rewind();
  return eventList;
}


void PwaGen::generate(std::shared_ptr<AbsLh> theLh, fitParCol& theFitParams){

  Info << "\n******** PAWIAN Event Generator **********\n" << endmsg;

  EvtMTRandomEngine myRandom(GlobalEnv::instance()->parser()->randomSeed());
  EvtRandom::setRandomEngine(&myRandom);
  bool generateEvents=true;
  int noOfAcceptedEvts=0;
  int noOfEvtsToGenerate=GlobalEnv::instance()->parser()->noOfGenEvts();
  int noOfIterations=0;
  int currentEvtNo=0;
  _maxFitWeight=0.;

  theLh->updateFitParams(theFitParams);

  while(generateEvents){
    noOfIterations++;
    
    if(((!_genWithModel || _useEvtWeight) && noOfIterations == 1) || 
       (_genWithModel && !_useEvtWeight && noOfIterations == 2)){
       Info << "Starting event generation" << endmsg;
    }
    if(_genWithModel && !_useEvtWeight && noOfIterations == 1){
       Info << "Getting max weight" << endmsg;
    }   

    Event* anEvent;
    std::shared_ptr<EventList> currentEvtList = GeneratePspEventList(EVENTS_PER_ITERATION);
    std::shared_ptr<EvtDataBaseList> eventListPtr(new EvtDataBaseList(0));
    
    while ((anEvent = currentEvtList->nextEvent()) !=0){
      currentEvtNo++;
      std::shared_ptr<EvtData> currentEvtData(eventListPtr->convertEvent(anEvent, currentEvtNo));

      if(!_genWithModel){
	  DumpEventToAsciiFile(currentEvtData); 
	  noOfAcceptedEvts++;
      }
      else if(_genWithModel && _useEvtWeight){
	 double fitWeight= theLh->calcEvtIntensity(currentEvtData.get(), theFitParams);
	 DumpEventToAsciiFile(currentEvtData, fitWeight);
	 noOfAcceptedEvts++;
      }
      else if(_genWithModel && !_useEvtWeight){
	 double fitWeight= theLh->calcEvtIntensity(currentEvtData.get(), theFitParams);
	 UpdateMaxFitWeight(fitWeight, noOfIterations);
	 
	 if(noOfIterations > 1){
	    double randWeight = EvtRandom::Flat(0., _maxFitWeight );
	    if (randWeight < fitWeight){
	       DumpEventToAsciiFile(currentEvtData);
	       noOfAcceptedEvts++;
	    }
	 }	 
      }

      if(noOfAcceptedEvts==noOfEvtsToGenerate){
	 generateEvents=false;
	 break;
      }
    }

    Info << "Iteration " << noOfIterations << " finished. Accepted events:\t" << noOfAcceptedEvts << endmsg;

    currentEvtList->rewind();
    currentEvtList->removeAndDeleteEvents(0,currentEvtList->size()); 
  }
}



void PwaGen::AddEventToEventList(std::shared_ptr<EventList> evtList, EvtVector4R* evt4Vec4Rs,int evtNumber, double weight){
  Event* newEvent = new Event();
  newEvent->addWeight(weight);
  for(unsigned int t=0; t<_finalStateParticles.size(); ++t){
    newEvent->addParticle(evt4Vec4Rs[t].get(0), evt4Vec4Rs[t].get(1), evt4Vec4Rs[t].get(2), evt4Vec4Rs[t].get(3));
  }

  evtList->add(newEvent);
}


void PwaGen::DumpEventToAsciiFile(std::shared_ptr<EvtData> evtData, double weight){

  std::vector<Particle* > fsParticles = GlobalEnv::instance()->Channel()->finalStateParticles();
  std::vector<Vector4<double>> current4Vecs;

  if(_useEvtWeight && _genWithModel){
      *_stream << weight << std::endl;
  }

  for(auto fspIter = fsParticles.begin(); fspIter != fsParticles.end(); ++fspIter ) {
    Vector4<double> tmp4vec = evtData->FourVecsString[ (*fspIter)->name() ];
    current4Vecs.push_back(tmp4vec);
    if(_energyFirst){
      *_stream << std::setprecision(8)  << tmp4vec.E()*_unitScaleFactor  << "\t" 
	       << tmp4vec.Px()*_unitScaleFactor << "\t" << tmp4vec.Py()*_unitScaleFactor << "\t" 
	       << tmp4vec.Pz()*_unitScaleFactor << "\t" << std::endl;
    }
    else{
      *_stream << std::setprecision(8)  << tmp4vec.Px()*_unitScaleFactor << "\t" 
	       << tmp4vec.Py()*_unitScaleFactor << "\t" << tmp4vec.Pz()*_unitScaleFactor << "\t" 
	       << tmp4vec.E()*_unitScaleFactor << std::endl;
    }
  }

  inv01MassWeightH1->Fill((current4Vecs.at(0)+current4Vecs.at(1)).M(), weight);
  if(_finalStateParticles.size()>2){
    inv02MassWeightH1->Fill((current4Vecs.at(0)+current4Vecs.at(2)).M(), weight);
    inv12MassWeightH1->Fill((current4Vecs.at(1)+current4Vecs.at(2)).M(), weight);
  }
}



void PwaGen::UpdateMaxFitWeight(double weight, int currentIteration){
   if(_maxFitWeight < weight){
      _maxFitWeight = weight;

      if(currentIteration <= 1){
	 Info << "Current max weight = " << _maxFitWeight << endmsg;
      }
      else{
	 Warning << "Raised max weight to " << _maxFitWeight << " while generating!" << endmsg;
      }
   }
}
