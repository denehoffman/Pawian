//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// PiPiScatteringXdecAmps class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/PiPiScatteringXdecAmps.hh"
#include "PwaUtils/PiPiScatteringDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/DynRegistry.hh"
#include "PwaUtils/TMatrixDynamics.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"
#include "FitParams/AbsPawianParameters.hh"

PiPiScatteringXdecAmps::PiPiScatteringXdecAmps(std::shared_ptr<PiPiScatteringDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
{
  initialize();
}

PiPiScatteringXdecAmps::PiPiScatteringXdecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
 AbsXdecAmp(theDec, channelID)
{
}

PiPiScatteringXdecAmps::~PiPiScatteringXdecAmps()
{
}

complex<double> PiPiScatteringXdecAmps::XdecPartAmp(const Spin& lamX, Spin& lamDec, short fixDaughterNr,
						    EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp) {
  complex<double> result(0.,0.);
  return result;
}

complex<double> PiPiScatteringXdecAmps::XdecAmp(const Spin& lamX, EvtData* theData, AbsXdecAmp* grandmaAmp) {
  complex<double> result(0.,0.);
  _tMatDyn->eval(theData, grandmaAmp, _orbMom);
  return result;
}

void PiPiScatteringXdecAmps::doTcheck(){
  _tMatDyn->doTcheck(_orbMom);
}

void PiPiScatteringXdecAmps::calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp) {
  AbsXdecAmp::calcDynamics(theData, grandmaAmp);
 return;
}

void  PiPiScatteringXdecAmps::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar) {
  _absDyn->fillDefaultParams(fitPar);
}

void PiPiScatteringXdecAmps::fillParamNameList() {
  _absDyn->fillParamNameList();
}

void PiPiScatteringXdecAmps::print(std::ostream& os) const {
  return; //dummy
}


void PiPiScatteringXdecAmps::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar) {
  _absDyn->updateFitParams(fitPar);
}

void PiPiScatteringXdecAmps::initialize() {
  _absDyn = DynRegistry::instance()->getDynamics(_decay);
  if(!(_decay->dynType() == "FVectorCompare"))
    _absDyn->addGrandMa(_decay);
  _tMatDyn=std::dynamic_pointer_cast<TMatrixDynamics>(_absDyn);
  _orbMom=_tMatDyn->orbitalL();
}
