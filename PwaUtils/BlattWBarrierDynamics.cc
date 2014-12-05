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
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/BlattWBarrierDynamics.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/BreitWignerFunction.hh"
#include "ConfigParser/ParserBase.hh"

BlattWBarrierDynamics::BlattWBarrierDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, const std::string& wignerDKey, double qR) :
  AbsDynamics(name, fsParticles, mother)
  ,_wignerDKey(wignerDKey)
  ,_qR(qR)
  ,_fitqRVals(false)
  ,_fitqRKey(_massKey+"qRPosOther")
{
  Info << "BlattWBarrierDynamics for " << _name <<endmsg;
  if(GlobalEnv::instance()->parser()->fitqRProduction()) _fitqRVals=true;
}

BlattWBarrierDynamics::~BlattWBarrierDynamics()
{
}

complex<double> BlattWBarrierDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  complex<double> result(1.,0.);
  if(OrbMom==0) return result;
  result=BarrierFactor::BlattWeisskopf(OrbMom, theData->DoubleString.at(_wignerDKey), _qR) /
    BarrierFactor::BlattWeisskopf(OrbMom, theData->DoubleString.at(_wignerDKey + "qNorm"), _qR);

  return result;
}

void  BlattWBarrierDynamics::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
  if(!_fitqRVals) return;
  fitVal.otherParams[_fitqRKey]=_qR;
  fitErr.otherParams[_fitqRKey]=0.3;
}

bool BlattWBarrierDynamics::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;
 
  if(_fitqRVals){
    if (!CheckDoubleEquality( _qR, theParamVal.otherParams.at(_fitqRKey))) _recalculate=true;
  }

  return _recalculate;
}

void BlattWBarrierDynamics::updateFitParams(fitParams& theParamVal){
  if(!_fitqRVals) return;
  _qR=theParamVal.otherParams.at(_fitqRKey);
}

void BlattWBarrierDynamics::setMassKey(std::string& theMassKey){
  _massKey=theMassKey;
  _fitqRKey=_massKey+"qRPosOther";
}


