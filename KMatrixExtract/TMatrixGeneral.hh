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
//TMatrixGeneral class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <memory>

#include "TROOT.h"


class TFile;
class TH1F;
class TH2F;
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

class TMatrixGeneral {

public:

  // create/copy/destroy:

  ///Constructor 
  TMatrixGeneral(pipiScatteringParser* theParser);


  /** Destructor */
  virtual ~TMatrixGeneral();

  virtual void fillParams();
  virtual void process();
  virtual void initHistos();

  // Getters:
 
protected:
  pipiScatteringParser* _pipiScatteringParser;
  std::shared_ptr<PiPiScatteringChannelEnv> _pipiScatteringChannelEnv; 
  std::string _projectionParticleNames;
  std::string _motherParticleName;
  int _decProjectionIndex;
  unsigned int _noOfSteps;
  double _stepSize; 
  double _massMin;
  double _massMax;
  std::vector< std::string> _gFactorNames;
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
  int _numStepsForSheetScan;
  Particle* _motherParticle;
  std::vector<Particle*> _fsParticles;
  TFile* _theTFile;

private:
  std::vector<TH1F*> _AmpRealH1Vec;
  std::vector<TH1F*> _AmpImagH1Vec;
  std::vector<TH1F*> _ImagT11m1H1Vec;
  std::vector<TH2F*> _ArgandH2Vec;
  std::vector<TH2F*> _PhaseH2Vec;
  std::vector<TH1F*> _ElasticityH1Vec;
  std::vector<TH1F*> _SqrT11H1Vec;
  std::vector<TH1F*> _phpH1Vec;
  std::vector<TH1F*> _phpH1RealVec;
  std::vector<TH1F*> _phpH1ImagVec;
  std::vector<TH1F*> _SqrtT1iH1Vec;
  std::vector<TH1F*> _ArgandUnits1iH1Vec;
  std::vector<TH1F*> _delta1iVec;
  std::vector<TH1F*> _SqrS1iH1Vec;
  std::vector<TH1F*> _speedPlotH1Vec;

  void init();
};


