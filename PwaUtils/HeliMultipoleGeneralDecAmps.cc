//************************************************************************//
//									  //
//  Copyright 2022 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// HeliMultipoleGeneralDecAmps class definition file. -*- C++ -*-
// Copyright 2022 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <sstream>

#include "PwaUtils/HeliMultipoleGeneralDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"
#include "FitParams/AbsPawianParameters.hh"

HeliMultipoleGeneralDecAmps::HeliMultipoleGeneralDecAmps(std::shared_ptr<IsobarHeliDecay> theDec, ChannelID channelID) :
  HeliMultipoleDecNonRefAmps(theDec, channelID)
{
}

HeliMultipoleGeneralDecAmps::HeliMultipoleGeneralDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  HeliMultipoleDecNonRefAmps(theDec, channelID)
{
}

HeliMultipoleGeneralDecAmps::~HeliMultipoleGeneralDecAmps()
{
}

complex<double> HeliMultipoleGeneralDecAmps::heliAmpLoop(EvtData* theData, const Spin& lamX, Spin& lam1, Spin& lam2, const Spin& J, bool isSym){
  complex<double> result(0.,0.);
  Spin lambda = lam2-lam1;
  Spin lam1Sym = std::abs(lam1);
  Spin lam2Sym = std::abs(lam2);
  Spin lambdaSym = lam2Sym-lam1Sym;
  bool doCalc=true;
  if( fabs(lambda) > J) doCalc=false;
  
  if(_daughter1IsStable && (_lam1MinProj!=lam1)) doCalc=false;
  if(_daughter2IsStable && (_lam2MinProj!=lam2)) doCalc=false;
  // InfoMsg << "J: " << J << " lamX: " << lamX << " lam1: " << lam1 << " lam2: " << lam2 << " lambda: " << lambda << " doCalc: " << doCalc << endmsg;
  if(doCalc){
   unsigned int IdJLamXLam12=FunctionUtils::spin3Index(_J, lamX, lambda);
   complex<double> currentAmp(0.,0.);
   double currentParity=1.;
   if(isSym) currentParity=(-1.)*_daughter2->theParity()*_JPCPtr->P*pow(-1,(1+_daughter2->J()-_JPCPtr->J));
   for (int i=0; i<_noOfAmps; ++i){
     //if(int(_JgammaMap.at(i))%2 == 0 && lam1<0) parityFactor *= -1.; //this is a fix which is needed to be checked
     currentAmp+=sqrt((2.*_JgammaMap.at(i)+1.)/(2.*_daughter2->J()+1.))
       *currentParity
       *Clebsch(_JgammaMap.at(i), lam1Sym,_JPCPtr->J, lambdaSym,  _daughter2->J(), lam2Sym)
       *_currentParamLocalMagExpi.at(i);
   }
   result+=currentAmp*theData->WignerDIdId3.at(_decay->wigDWigDRefId()).at(IdJLamXLam12) *daughterAmp(lam1, lam2, theData);
   // InfoMsg<< "currentAmp: " << currentAmp << endmsg;
   // InfoMsg<< "theData->WignerDIdId3.at(_decay->wigDWigDRefId()).at(IdJLamXLam12): " << theData->WignerDIdId3.at(_decay->wigDWigDRefId()).at(IdJLamXLam12) << endmsg;
   // InfoMsg<< "result: " << result << endmsg;
  }
  return result;
}
