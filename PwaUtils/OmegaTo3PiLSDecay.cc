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

// OmegaTo3PiLSDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "PwaUtils/OmegaTo3PiLSDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "Utils/FunctionUtils.hh"

OmegaTo3PiLSDecay::OmegaTo3PiLSDecay(Particle* mother, Particle* daughter1, Particle* daughter2, Particle* daughter3, ChannelID channelID) :
  OmegaTo3PiDecay(mother, daughter1, daughter2, daughter3, channelID)
{
  _wignerDKey="omegaTo3piLS_"+FunctionUtils::particleListName(_finalStateParticles);
}

OmegaTo3PiLSDecay::~OmegaTo3PiLSDecay(){
}

void OmegaTo3PiLSDecay::fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData){
  Vector4<double> all4Vec(0.,0.,0.,0.);
  Vector4<double> mother4Vec(0.,0.,0.,0.);

  //fill all4Vec
  std::map<std::string, Vector4<double> >::iterator itMap;
  for(itMap=fsMap.begin(); itMap!=fsMap.end(); ++itMap){
    all4Vec+=itMap->second;
  }

  //fill mother4Vec
  std::vector<Particle*>::iterator itP;
  for(itP = _finalStateParticles.begin(); itP != _finalStateParticles.end(); ++itP){
    itMap=fsMap.find((*itP)->name());
    mother4Vec+=itMap->second;
  }

  //fill daughter1 and daughter2 4Vec
  itMap=fsMap.find(_daughter1->name());
  Vector4<double> daughter1_4Vec=itMap->second;
  itMap=fsMap.find(_daughter2->name());
  Vector4<double> daughter2_4Vec=itMap->second;
  itMap=fsMap.find(_daughter3->name());
  Vector4<double> daughter3_4Vec=itMap->second;

  //  Vector4<double> daughter1_HeliOmega=helicityVec(all4Vec, mother4Vec, daughter1_4Vec);
  //  Vector4<double> daughter2_HeliOmega=helicityVec(all4Vec, mother4Vec, daughter2_4Vec);
  //  Vector4<double> daughter3_HeliOmega=helicityVec(all4Vec, mother4Vec, daughter3_4Vec);

  Vector4<double> daughter1_HeliOmega=helicityVec(prodParticle4Vec, mother4Vec, daughter1_4Vec);
  Vector4<double> daughter2_HeliOmega=helicityVec(prodParticle4Vec, mother4Vec, daughter2_4Vec);
  Vector4<double> daughter3_HeliOmega=helicityVec(prodParticle4Vec, mother4Vec, daughter3_4Vec);

  //calculate normal of the decay plane
  Vector4<float> normOmegaDecHeliOmega_4V(0.5*(daughter1_HeliOmega.T()+daughter2_HeliOmega.T()+daughter3_HeliOmega.T()),
					  daughter1_HeliOmega.Y()*daughter2_HeliOmega.Z()-daughter1_HeliOmega.Z()*daughter2_HeliOmega.Y(),
					  daughter1_HeliOmega.Z()*daughter2_HeliOmega.X()-daughter1_HeliOmega.X()*daughter2_HeliOmega.Z(),
					  daughter1_HeliOmega.X()*daughter2_HeliOmega.Y()-daughter1_HeliOmega.Y()*daughter2_HeliOmega.X());

  double theQ=daughter1_HeliOmega.E()-daughter1_HeliOmega.M()+daughter2_HeliOmega.E()-daughter2_HeliOmega.M()+daughter3_HeliOmega.E()-daughter3_HeliOmega.M();
  double lambdaNorm=theQ*theQ*(theQ*theQ/108.+daughter1_HeliOmega.M()*theQ/9.+daughter1_HeliOmega.M()*daughter1_HeliOmega.M()/3.);
  double lambdaOmegaDec=normOmegaDecHeliOmega_4V.P()*normOmegaDecHeliOmega_4V.P()/lambdaNorm;
  evtData->DoubleId[_wigDWigDRefId]=lambdaOmegaDec;

  for(Spin lamOmega=-1;  lamOmega<=1; ++lamOmega){
    Id1StringType IdLamOmega=FunctionUtils::spin1Index(lamOmega);
     evtData->WignerDIdId1[_wigDWigDRefId][IdLamOmega]=Wigner_D(normOmegaDecHeliOmega_4V.Phi(), normOmegaDecHeliOmega_4V.Theta(),0, 1,lamOmega,0);
  }

}


