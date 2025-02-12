//************************************************************************//
//									  //
//  Copyright 2025 Bertram Kopf (bertram@ep1.rub.de)			  //
//          	   - Ruhr-Universit??t Bochum 				  //
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

// SigmaParamDynamics class definition file. -*- C++ -*-
// Copyright 2025 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/SigmaParamDynamics.hh"
#include "Utils/PawianIOUtils.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "PwaDynamics/SigmaParameterization.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "Utils/IdStringMapRegistry.hh"
#include "Utils/PawianConstants.hh"

SigmaParamDynamics::SigmaParamDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  TMatrixSigmaParDynamics(name, fsParticles, mother)
   ,_currentStrength(std::complex<double>(1.,0.))  
  ,_strenghtNameMag(_massKey+"StrengthMag")
  ,_strenghtNamePhi(_massKey+"StrengthPhi")
{
}

SigmaParamDynamics::~SigmaParamDynamics()
{
}

complex<double> SigmaParamDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

  if ( !_recalculate){
    return _cachedMap[theData->evtNo];
  }

  complex<double> result=_sigmaParameterization->calc(theData->DoubleMassId.at(_dynId), _currentStrength);
  if ( _cacheAmps){
     _cachedMap[theData->evtNo]=result;
  }  

  return result;
}

void SigmaParamDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  TMatrixSigmaParDynamics::fillDefaultParams(fitPar);
  //fill production strength mag
  fitPar->Add(_strenghtNameMag, std::abs(_currentStrength), 0.01);
  //fill production strength phase
  fitPar->Add(_strenghtNamePhi, std::arg(_currentStrength), 0.01);
}

void SigmaParamDynamics::fillParamNameList(){
  _paramNameList.clear();
  TMatrixSigmaParDynamics::fillParamNameList();
  _paramNameList.push_back(_strenghtNameMag);
  _paramNameList.push_back(_strenghtNamePhi);
}

bool SigmaParamDynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){
  return AbsParamHandler::checkRecalculation(fitParNew, fitParOld);
}

void SigmaParamDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){

  TMatrixSigmaParDynamics::updateFitParams(fitPar);
  double currentStrengthMag = fitPar->Value(_strenghtNameMag);
  double currentStrengthPhi = fitPar->Value(_strenghtNamePhi);
  _currentStrength=std::polar(currentStrengthMag,currentStrengthPhi);
}



