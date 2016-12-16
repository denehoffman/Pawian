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
#include "Utils/IdStringMapRegistry.hh"

BlattWBarrierDynamics::BlattWBarrierDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, const std::string& wignerDKey, double qR) :
  AbsDynamics(name, fsParticles, mother)
  ,_wignerDKey(wignerDKey)
  ,_wignerDqNormKey(_wignerDKey+"qNorm")
  ,_wignerDqId(IdStringMapRegistry::instance()->keyStringId(_keyForMassList, _wignerDKey))
  ,_wignerDqNormId(IdStringMapRegistry::instance()->keyStringId(_keyForMassList, _wignerDqNormKey))
  ,_qR(qR)
  ,_fitqRVals(false)
  ,_fitqRKey(_massKey+"qRPosOther")
{
  InfoMsg << "BlattWBarrierDynamics for " << _name <<endmsg;
  if(GlobalEnv::instance()->parser()->fitqRProduction()) _fitqRVals=true;
  _isLdependent=true;
  // _qR should be between 0.02  and 20. 
  if (_qR<0.02){
    WarningMsg << "_qR value of " << _qR << " to low! Set it to 0.02!!!" << endmsg;
    _qR=0.02;
  }
  if (_qR>20.){
    WarningMsg << "_qR value of " << _qR << " to high! Set it to 20.!!!" << endmsg;
    _qR=20.;
  }
}

BlattWBarrierDynamics::~BlattWBarrierDynamics()
{
}

complex<double> BlattWBarrierDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  complex<double> result(1.,0.);
  if(OrbMom==0) return result;
  result=BarrierFactor::BlattWeisskopf(OrbMom, theData->DoubleMassId.at(_wignerDqId), _qR) /
    BarrierFactor::BlattWeisskopf(OrbMom, theData->DoubleMassId.at(_wignerDqNormId), _qR);

  return result;
}

void  BlattWBarrierDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  if(!_fitqRVals) return;
  fitPar->Add(_fitqRKey, _qR, 0.05);
  fitPar->SetLimits(_fitqRKey, 0.02, 20.);
}

void BlattWBarrierDynamics::fillParamNameList(){
  if(!_fitqRVals) return;
  _paramNameList.clear();
  _paramNameList.push_back(_fitqRKey);
}

void BlattWBarrierDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  if(!_fitqRVals) return;
  _qR = fitPar->Value(_fitqRKey);
}

void BlattWBarrierDynamics::setMassKey(std::string& theMassKey){
  _massKey=theMassKey;
  _fitqRKey=_massKey+"qRPosOther";
}


