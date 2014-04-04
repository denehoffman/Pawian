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

BreitWignerRelDynamics::BreitWignerRelDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, double massSumDaughter1, double massSumDaughter2) :
  BreitWignerDynamics(name, fsParticles, mother)
{
  _fsp1Mass=massSumDaughter1;
  _fsp2Mass=massSumDaughter2;
  Info << "BreitWignerRelDynamics for " << _name << " with  massSumDaughter1= " << _fsp1Mass << " and  massSumDaughter2= " << _fsp2Mass <<endmsg;
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

