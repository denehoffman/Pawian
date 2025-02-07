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

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <memory>

#include "PwaUtils/AbsDynamics.hh"
#include "FitParams/AbsPawianParameters.hh"

class SigmaParameterization;

class SigmaParamDynamics : public AbsDynamics{

public:
  SigmaParamDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother);
  virtual ~SigmaParamDynamics();

  virtual std::string type() {return "SigmaParamDynamics";}
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom=0);
  
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);

  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();
protected:
  //  std::string _massKey;
  std::shared_ptr<SigmaParameterization> _sigmaParamDyn; 
  complex<double> _currentStrength;
  double _currentPoleMass;
  double _currentb1;
  double _currentb2;
  double _currenta;
  double _currentg4pi;
  const std::string _strenghtNameMag;
  const std::string _strenghtNamePhi;
  const std::string _poleMassName;
  const std::string  _b1Name;
  const std::string  _b2Name;
  const std::string  _aName;
  const std::string  _g4piName;
  
private:

};
