//************************************************************************//
//									  //
//  Copyright 2021 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// VoigtBlattWRelDynamics class definition file. -*- C++ -*-
// Copyright 2021 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <memory>

#include "PwaUtils/BreitWignerBlattWRelDynamics.hh"
#include "PwaDynamics/BarrierFactor.hh"

class VoigtBlattWRelDynamics : public BreitWignerBlattWRelDynamics{

public:
  VoigtBlattWRelDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::vector<Particle*>& fsParticlesDaughter1, std::vector<Particle*>& fsParticlesDaughter2, unsigned int prodOrbMom, double qR=BarrierFactor::qRDefault);
  virtual ~VoigtBlattWRelDynamics();

  virtual std::string type() {return "VoigtBlattWRelDynamics";}
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom=0);
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();   

  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void setMassKey(std::string& theMassKey);
  
protected:
  std::string _massSigmaKey;
  double _currentSigma;
  int _prodOrbMom;
  
private:

};
