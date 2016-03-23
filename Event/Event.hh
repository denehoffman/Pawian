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

#pragma once

#include "qft++/topincludes/tensor.hh"
#include "qft++/tensor/Vector4.hh"
#include "ErrLogger/ErrLogger.hh"
#include <vector>

enum pidType {proton, pion, kaon, muon, electron, photon, resonance, maxPid};

struct EvtPartData {
  Vector4<float>* vector4;
  float pidVector[maxPid];
  EvtPartData(double e, double px, double py, double pz) {
    vector4 = new Vector4<float>(e,px,py,pz);
   }
  EvtPartData(const EvtPartData& evtPartData) {
    vector4 = new Vector4<float>(evtPartData.vector4);
    for(unsigned int i=0; i<maxPid; ++i) pidVector[i]=evtPartData.pidVector[i];
   }
  ~EvtPartData() {
    delete vector4;
  }
};

class Event
{
public:
  Event();
  Event(int evtNumber);
  ~Event();

  void addParticle(double e, double px, double py, double pz);
  void addWeight(double);

  Vector4<float>* p4(unsigned int i);
  float* pid(unsigned int i);

  int size();
  double Weight(){ return weight; }
  int eventNo() const {return evtNo; }
  
  bool operator<(const Event& compare) const;

private:
  std::vector<EvtPartData*> particles;
  double weight;
  int evtNo;
};
