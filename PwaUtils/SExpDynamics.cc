//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// SExpDynamics class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/SExpDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

SExpDynamics::SExpDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& daughterName) :
  AbsDynamics(name, fsParticles, mother)
  ,_currentL(0.)
  ,_daughterName(daughterName)
{
  _isLdependent=false;
}

SExpDynamics::~SExpDynamics()
{
}

complex<double> SExpDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  //InfoMsg << "theData->DoubleMassId.at(_dynId): " << theData->DoubleMassId.at(_dynId) << endmsg;
  complex<double> result(pow(theData->DoubleMassId.at(_dynId), _currentL),0.);
  
  return result;
}


void SExpDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  //fill mass
  std::string massName=_massKey+_daughterName+"L";
  
  fitPar->Add(massName, 0., 0.01);
  fitPar->SetLimits(massName, -3., 3.);
}


void SExpDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  std::string massName=_massKey+_daughterName+"L";
  _currentL=fitPar->Value(massName);
}

void SExpDynamics::fillParamNameList(){
  _paramNameList.clear();
  //fill mass
  std::string massName=_massKey+_daughterName+"L";
  _paramNameList.push_back(massName);
}
