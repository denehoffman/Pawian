//************************************************************************//
//									  //
//  Copyright 2023 Bertram Kopf (bertram@ep1.rub.de)			  //
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
//TMatrixGeneralBase class definition file. -*- C++ -*-
// Copyright 2023 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <memory>

class KMatrixParser;
class AbsPhaseSpace;
class TMatrixDynamics;
class TMatrixRel;
class KMatrixRel;
class KPole;
class ParticleTable;
class Particle;
class KMatrixParser;
class AbsPawianParameters;
class pipiScatteringParser;
class PiPiScatteringChannelEnv;

class TMatrixGeneralBase {

public:

  // create/copy/destroy:

  ///Constructor
  TMatrixGeneralBase(); 
  TMatrixGeneralBase(pipiScatteringParser* theParser);


  /** Destructor */
  virtual ~TMatrixGeneralBase();

  virtual void fillParams();
  std::vector<std::shared_ptr<AbsPhaseSpace> > getPhps() {return _phpVecs;}
  virtual void SetParamValue(const std::string & paramName, double paramVal);
  virtual double GetParamValue(const std::string & paramName);
  // Getters:

protected:
  pipiScatteringParser* _pipiScatteringParser;
  std::shared_ptr<PiPiScatteringChannelEnv> _pipiScatteringChannelEnv; 
  std::string _projectionParticleNames;
  std::string _motherParticleName;
  int _decProjectionIndex;

  std::vector< std::string> _gFactorNames;
  std::vector<std::string> _kMatrixParamNames;
  std::vector<std::shared_ptr<AbsPhaseSpace> > _phpVecs;
  std::shared_ptr<AbsPhaseSpace> _phpVecCurrent;
  ParticleTable* _particleTable;
  std::shared_ptr<AbsPawianParameters> _params;
  std::string _pathToFitParams;
  int _orbitalL;
  std::string _pathToKMatrixParser;
  std::shared_ptr<KMatrixParser> _kMatrixParser;
  std::shared_ptr<TMatrixDynamics> _tMatrDyn;
  std::shared_ptr<TMatrixRel> _tMatr;
  std::shared_ptr<KMatrixRel> _kMatr;
  std::vector<double> _energyPlaneBorders;
  Particle* _motherParticle;
  std::vector<Particle*> _fsParticles;
  double _massMin;
  double _massMax;
private:
  void init();
};


