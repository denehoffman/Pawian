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

// BreitWignerDynamics class definition file. -*- C++ -*-
// Copyright 20123Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/BreitWignerRelDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/BreitWignerFunction.hh"

BreitWignerRelDynamics::BreitWignerRelDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  BreitWignerDynamics(name, fsParticles, mother)
{
  Particle* fsParticle1=fsParticles[0];
  Particle* fsParticle2=fsParticles[1];
  if(0==fsParticle1 || 0==fsParticle2){
    Alert << "0 pionter for final state particle 0 or 1!!!" << endmsg;
    exit(0);
  }
  _fsp1Mass=fsParticle1->mass();
  _fsp2Mass=fsParticle2->mass();
}

BreitWignerRelDynamics::~BreitWignerRelDynamics()
{
}

complex<double> BreitWignerRelDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate){
    return _cachedMap[evtNo];
  }

  complex<double> result=BreitWignerFunction::Rel(theData->DoubleString.at(_dynKey), _currentMass, _currentWidth, _fsp1Mass, _fsp2Mass);  
  if ( _cacheAmps){
     theMutex.lock();
     _cachedMap[evtNo]=result;
     theMutex.unlock();
  }  

  return result;
}

