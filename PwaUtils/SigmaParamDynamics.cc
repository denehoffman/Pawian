//************************************************************************//
//									  //
//  Copyright 2025 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// SigmaParamDynamics class definition file. -*- C++ -*-
// Copyright 2025 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/SigmaParamDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/SigmaParameterization.hh"

SigmaParamDynamics::SigmaParamDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsDynamics(name, fsParticles, mother)
  ,_sigmaParamDyn(std::shared_ptr<SigmaParameterization>(new SigmaParameterization))
  ,_currentStrength(std::complex<double>(1.,0.))  
  ,_currentPoleMass(_sigmaParamDyn->barePoleMass())
  ,_currentb1(_sigmaParamDyn->b1Param())
  ,_currentb2(_sigmaParamDyn->b2Param())
  ,_currenta(_sigmaParamDyn->aParam())
  ,_currentg4pi(_sigmaParamDyn->g4piParam())
  ,_strenghtNameMag(_massKey+"StrengthMag")
  ,_strenghtNamePhi(_massKey+"StrengthPhi")
  ,_poleMassName("sigmaPoleBareMass")
  ,_b1Name("sigmab1")
  ,_b2Name("sigmab2")
  ,_aName("sigmaa")
  ,_g4piName("sigmag4pi")
{
  _isLdependent=false;
}

SigmaParamDynamics::~SigmaParamDynamics()
{
}

complex<double> SigmaParamDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  if ( !_recalculate){
    return _cachedMap[theData->evtNo];
  }

  complex<double> result=_sigmaParamDyn->calc(theData->DoubleMassId.at(_dynId), _currentStrength);
  if ( _cacheAmps){
     _cachedMap[theData->evtNo]=result;
  }  

  return result;
}


void SigmaParamDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  //fill production strength mag
  fitPar->Add(_strenghtNameMag, std::abs(_currentStrength), 0.01);
  //fill production strength phase
  fitPar->Add(_strenghtNamePhi, std::arg(_currentStrength), 0.01);
  //fill pole bare mass
  fitPar->Add(_poleMassName, _currentPoleMass, 0.01, std::abs(_currentPoleMass-0.2), _currentPoleMass+0.2);
  //fill b1 param
  fitPar->Add(_b1Name, _currentb1, 0.01, std::abs(_currentb1-0.2), _currentb1+0.2);
  //fill b2 param
  fitPar->Add(_b2Name, _currentb2, 0.01, std::abs(_currentb2-0.2), _currentb2+0.2);
  //fill a param
  fitPar->Add(_aName, _currenta, 0.01, std::abs(_currenta-0.2), _currenta+0.2);
  //fill g4pi param
  fitPar->Add(_g4piName, _currentg4pi, 0.00001, 0., _currentg4pi+0.0001);
}


void SigmaParamDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  double currentStrengthMag = fitPar->Value(_strenghtNameMag);
  double currentStrengthPhi = fitPar->Value(_strenghtNamePhi);
  _currentStrength=std::polar(currentStrengthMag,currentStrengthPhi);

  _currentPoleMass = fitPar->Value(_poleMassName);
  _sigmaParamDyn->setBarePoleMass(_currentPoleMass);

  _currentb1= fitPar->Value(_b1Name);
  _sigmaParamDyn->setb1Param(_currentb1);

  _currentb2= fitPar->Value(_b2Name);
  _sigmaParamDyn->setb2Param(_currentb2);

  _currenta= fitPar->Value(_aName);
  _sigmaParamDyn->setaParam(_currenta);

  _currentg4pi= fitPar->Value(_g4piName);
  _sigmaParamDyn->setg4piParam(_currentg4pi);  
}

void SigmaParamDynamics::fillParamNameList(){
  _paramNameList.clear();
  //fill
  _paramNameList.push_back(_strenghtNameMag);
  _paramNameList.push_back(_strenghtNamePhi);
  _paramNameList.push_back(_poleMassName);
  _paramNameList.push_back(_b1Name);
  _paramNameList.push_back(_b2Name);
  _paramNameList.push_back(_aName);
  _paramNameList.push_back(_g4piName);
}
