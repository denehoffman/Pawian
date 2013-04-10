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

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PspGen/EvtVector4R.hh"
#include "PwaUtils/FitParamsBase.hh"

#include "TROOT.h"

class AbsEnv;
class AbsLh;
class EvtData;
class EvtVector4R;
class Particle;
class TFile;
class TH1F;
class EventList;

class PwaGen{

public:
  PwaGen(AbsEnv* theEnv);
  virtual ~PwaGen();

  virtual void generate(boost::shared_ptr<AbsLh> theLh, fitParams& theFitParams);

protected:

  void addEvt(EventList& evtList, EvtVector4R* evt4Vec4R, int evtNumber, double weight=1.);

  //  void dumpAscii(EvtVector4R* evt4Vec4R);
  void dumpAscii(EvtData* evtData, double weight=1.);

  AbsEnv* _absEnv;
  EvtVector4R _initial4Vec;
  std::vector<Particle*> _finalStateParticles;
  
  std::ofstream* _stream;

  double mass[30];
  //  EvtVector4R p4[30];

  TFile* _theTFile;
  TH1F* inv01MassH1;
  TH1F* inv02MassH1;
  TH1F* inv12MassH1;
  bool _genWithModel;
  double _unitScaleFactor;
  bool _energyFirst;
  bool _useEvtWeight;
};
