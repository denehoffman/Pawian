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

// TMatrixSigmaParDynamics class definition file. -*- C++ -*-
// Copyright 2025 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <memory>

#include "PwaUtils/TMatrixDynamics.hh"

class AbsXdecAmp;
class AbsPhaseSpace;
class AbsPawianParameters;
class SigmaParameterization;

class TMatrixSigmaParDynamics : public TMatrixDynamics{

public:
  //  TMatrixSigmaParDynamics();
  TMatrixSigmaParDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser, std::string dataType, std::string projectionParticleNames="");
  TMatrixSigmaParDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother);
  virtual ~TMatrixSigmaParDynamics();

  virtual std::string type() {return "TMatrixSigmaParDynamics";}
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp=0, Spin OrbMom=0);
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();
  //  virtual bool checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld); 
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);

  int prodProjectionIndex() {return _prodProjectionIndex;}
  void setProdProjectionIndex(int idx);
  
protected:
  std::shared_ptr<SigmaParameterization> _sigmaParameterization;
  double _currentPoleMass;
  double _currentb1;
  double _currentb2;
  double _currenta;
  double _currentg4pi;
  //  const std::string _strenghtNameMag;
  //  const std::string _strenghtNamePhi;
  const std::string _poleMassName;
  const std::string  _b1Name;
  const std::string  _b2Name;
  const std::string  _aName;
  const std::string  _g4piName;  
  
  virtual void evalElasticity(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalPhase(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalRelativePhase(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalArgandUnits(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalTreal(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalTimag(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalTabs(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalPhasePhi(EvtData* theData, double currentMass, Spin OrbMom);

  
private:

};
