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
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <thread>

#include <boost/thread.hpp>

#include "PwaUtils/PwaGen.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "ConfigParser/ParserBase.hh"
#include "PwaUtils/GlobalEnv.hh"

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include "PspGen/EvtGenKine.hh"
#include "PspGen/EvtRandom.hh"
#include "PspGen/EvtMTRandomEngine.hh"
#include "PspGen/EvtVector4R.hh"

#include "Event/Event.hh"
#include "Event/EventList.hh"

#include "TFile.h"
#include "TH1F.h"

PwaGen::PwaGen() :
   _initial4Vec(EvtVector4R(GlobalEnv::instance()->Channel()->initial4Vec().E(), GlobalEnv::instance()->Channel()->initial4Vec().Px(),
                            GlobalEnv::instance()->Channel()->initial4Vec().Py(), GlobalEnv::instance()->Channel()->initial4Vec().Pz()))
  ,_finalStateParticles(GlobalEnv::instance()->Channel()->finalStateParticles())
  ,_genWithModel(true)
  ,_unitScaleFactor(1.)
  ,_energyFirst(false)
  ,_useEvtWeight(GlobalEnv::instance()->parser()->useEvtWeight())
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
  _genWithModel = GlobalEnv::instance()->parser()->generateWithModel();

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
}

PwaGen::~PwaGen()
{
  _theTFile->Write();
  _theTFile->Close();
  _stream->close();
}

void PwaGen::generate(std::shared_ptr<AbsLh> theLh, fitParams& theFitParams){
  // std::ofstream theStream ("output.dat");
  int counterFsp=0;
  std::vector<Particle*>::iterator it;
  for(it=_finalStateParticles.begin(); it!=_finalStateParticles.end(); ++it){
    mass[counterFsp]=(*it)->mass();
    counterFsp++;
  }

  EvtMTRandomEngine myRandom(GlobalEnv::instance()->parser()->randomSeed());
  EvtRandom::setRandomEngine(&myRandom);
  bool generateEvents=true;
  int noOfAcceptedEvts=0;
  int noOfEvtsToGenerate=GlobalEnv::instance()->parser()->noOfGenEvts();
  int noOfAllGenEvts=0;
  int noOfIterations=0;
  int currentEvtNo=0;
  double maxFitWeight=0.;

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
      if(_finalStateParticles.size()>2){
	inv02MassH1->Fill((p4[0]+p4[2]).mass());
	inv12MassH1->Fill((p4[1]+p4[2]).mass());
      }
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

    std::shared_ptr<EvtDataBaseList> eventListPtr(new EvtDataBaseList(0));

    std::vector<EvtData*> dataList;

    theLh->updateFitParams(theFitParams);
    while ((anEvent = currentEvtList.nextEvent()) !=0){
      EvtData* currentEvtData=eventListPtr->convertEvent(anEvent, currentEvtNo);
      currentEvtNo++;
  
      if(!_genWithModel){
	  dumpAscii(currentEvtData);
	  noOfAcceptedEvts++;
	  if(noOfAcceptedEvts==noOfEvtsToGenerate){
	    delete currentEvtData;
	    generateEvents=false;
	    break;
	  }
      }
      //    }
      else{ //generate with model
	//fit retrieve maxFitWeight
	theLh->updateFitParams(theFitParams);
	double fitWeight= theLh->calcEvtIntensity( currentEvtData, theFitParams );
	if (maxFitWeight< fitWeight) maxFitWeight=fitWeight;
	
	DebugMsg << currentEvtData->evtNo <<  "\tfitWeight:\t" << fitWeight << endmsg;
	if (_useEvtWeight){
	  dumpAscii(currentEvtData, fitWeight);
	  noOfAcceptedEvts++;
	  if(noOfAcceptedEvts==noOfEvtsToGenerate){
	    delete currentEvtData;
	    generateEvents=false;
	    break;
	  }
	}
	
	if (!_useEvtWeight && noOfIterations>1){
	  double randWeight = EvtRandom::Flat( 0., maxFitWeight );
	  if ( randWeight < fitWeight ){
	    dumpAscii(currentEvtData);
	    noOfAcceptedEvts++;
	    if(noOfAcceptedEvts==noOfEvtsToGenerate){
	      delete currentEvtData;
	      generateEvents=false;
	      break;
	    }
	  }
	}
      }
      delete currentEvtData;
      
    }
    Info << "iteration:\t" << noOfIterations << "\tnoOfAcceptedEvts:\t" << noOfAcceptedEvts << endmsg;

    //delete Events
    currentEvtList.rewind();
    currentEvtList.removeEvents(0,currentEvtList.size()); 
  }
  
}

void PwaGen::addEvt(EventList& evtList, EvtVector4R* evt4Vec4Rs,int evtNumber, double weight){
  Event* newEvent = new Event();
  newEvent->addWeight(weight);
  for(unsigned int t=0; t<_finalStateParticles.size(); ++t){
    newEvent->addParticle(evt4Vec4Rs[t].get(0), evt4Vec4Rs[t].get(1), evt4Vec4Rs[t].get(2), evt4Vec4Rs[t].get(3));
  }

  evtList.add(newEvent);
}


void PwaGen::dumpAscii(EvtData* evtData, double weight){
  if (_useEvtWeight && _genWithModel) *_stream << weight << std::endl;
  std::vector<Particle* > fsParticles = GlobalEnv::instance()->Channel()->finalStateParticles();

  std::vector<Vector4<double>> current4Vecs;

  std::vector<Particle* >::const_iterator fspIter = fsParticles.begin();
  for( ; fspIter != fsParticles.end(); ++fspIter ) {
    Vector4<double> tmp4vec = evtData->FourVecsString[ (*fspIter)->name() ];
    current4Vecs.push_back(tmp4vec);
    if(_energyFirst){
      *_stream << std::setprecision(8)  << tmp4vec.E()*_unitScaleFactor  << "\t" << tmp4vec.Px()*_unitScaleFactor << "\t" << tmp4vec.Py()*_unitScaleFactor << "\t" << tmp4vec.Pz()*_unitScaleFactor << "\t" << std::endl;
    }
    else{
      *_stream << std::setprecision(8)  << tmp4vec.Px()*_unitScaleFactor << "\t" << tmp4vec.Py()*_unitScaleFactor << "\t" << tmp4vec.Pz()*_unitScaleFactor << "\t" << tmp4vec.E()*_unitScaleFactor << std::endl;
    }

  }
  inv01MassWeightH1->Fill((current4Vecs[0]+current4Vecs[1]).M(), weight);
  if(_finalStateParticles.size()>2){
    inv02MassWeightH1->Fill((current4Vecs[0]+current4Vecs[2]).M(), weight);
    inv12MassWeightH1->Fill((current4Vecs[1]+current4Vecs[2]).M(), weight);
  }
}
