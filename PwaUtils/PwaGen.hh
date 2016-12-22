//************************************************************************//
//                                                                        //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)                      //
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
#include <memory>

#include "PspGen/EvtVector4R.hh"
//#include "FitParams/FitParColBase.hh"

#include "TROOT.h"

class AbsLh;
class EvtData;
class EvtVector4R;
class Particle;
class TFile;
class TH1F;
class EventList;
class AbsPawianParameters;




class PwaGen{

public:
  PwaGen();
  virtual ~PwaGen();
  virtual void generate(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> theFitParams);

protected:

  void AddEventToEventList(std::shared_ptr<EventList> evtList, EvtVector4R* evt4Vec4R, int evtNumber, double weight=1.);
  void DumpEventToAsciiFile(std::shared_ptr<EvtData> evtData, double weight=1.);
  void UpdateMaxFitWeight(double weight, int currentIteration);
  std::shared_ptr<EventList> GeneratePspEventList(unsigned int numEvents);

  bool _energyFirst;
  bool _useDataEvtWeight;
  bool _useMcEvtWeight;
  bool _genWithModel;
  const bool _useMassRange;
  const bool _usePhpDynamics;

  double _unitScaleFactor;
  double _fspMasses[30];
  double _maxFitWeight;

  EvtVector4R _initial4Vec;
  std::vector<Particle*> _finalStateParticles;
  std::ofstream* _stream;

  TFile* _theTFile;
  TH1F* inv01MassH1;
  TH1F* inv02MassH1;
  TH1F* inv12MassH1;
  TH1F* inv01MassWeightH1;
  TH1F* inv02MassWeightH1;
  TH1F* inv12MassWeightH1;
};
