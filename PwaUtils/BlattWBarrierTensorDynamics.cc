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

#include "PwaUtils/BlattWBarrierTensorDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

BlattWBarrierTensorDynamics::BlattWBarrierTensorDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, const std::string& wignerDKey, double qR) :
  BlattWBarrierDynamics(name, fsParticles, mother, wignerDKey, qR)
{
  _isLdependent=true;
}

BlattWBarrierTensorDynamics::~BlattWBarrierTensorDynamics()
{
}

complex<double> BlattWBarrierTensorDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  complex<double> result(1.,0.);
  if(OrbMom==0) return result;
  result=BarrierFactor::BlattWeisskopfTensor(OrbMom, theData->DoubleMassId.at(_wignerDqId), _qR) /
    BarrierFactor::BlattWeisskopfTensor(OrbMom, theData->DoubleMassId.at(_wignerDqNormId), _qR);

  return result;
}


