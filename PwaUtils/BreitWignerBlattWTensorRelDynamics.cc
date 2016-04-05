//************************************************************************//
//									  //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)			  //
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
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/BreitWignerBlattWTensorRelDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/BreitWignerFunction.hh"

BreitWignerBlattWTensorRelDynamics::BreitWignerBlattWTensorRelDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::vector<Particle*>& fsParticlesDaughter1, std::vector<Particle*>& fsParticlesDaughter2, double qR) :
  BreitWignerRelDynamics(name, fsParticles, mother, fsParticlesDaughter1, fsParticlesDaughter2)
  ,_qR(qR)
{
  _isLdependent=true;
}

BreitWignerBlattWTensorRelDynamics::~BreitWignerBlattWTensorRelDynamics()
{
  _isLdependent=true;
}

complex<double> BreitWignerBlattWTensorRelDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  int evtNo=theData->evtNo;
  int orbMom(OrbMom);
  if ( _cacheAmps && !_recalculate){
    return _cachedLMap.at(evtNo).at(orbMom);
  }

  complex<double> result=BreitWignerFunction::BlattWTensorRel(orbMom, theData->DoubleMassId.at(_dynId), _currentMass, _currentWidth, theData->DoubleMassId.at(_dynMassIdDaughter1), theData->DoubleMassId.at(_dynMassIdDaughter2), _qR);  
  if ( _cacheAmps){
     theMutex.lock();
     _cachedLMap[evtNo][orbMom]=result;
     theMutex.unlock();
  }  

  return result;
}


