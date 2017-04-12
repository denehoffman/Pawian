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
//MatrixPiPiSWaveSimple4piPhp class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "TROOT.h"


class TFile;
class TH1F;
class TH2F;
class AbsPhaseSpace;
class TMatrixRel;
class KMatrixBase;
class KPole;
class ParticleTable;

class MatrixPiPiSWaveSimple4piPhp {

public:

  // create/copy/destroy:

  ///Constructor 
  MatrixPiPiSWaveSimple4piPhp(int numStepsForSheetScan, std::vector<double> energyPlaneBorders);


  /** Destructor */
  virtual ~MatrixPiPiSWaveSimple4piPhp();

  // Getters:
 
protected:


private:
  const int _noOfSteps;
  double _stepSize;
  double _massMin;
  double _massMax;
  const int _noOfChannels;
   std::vector< std::string> _poleNames;
  std::vector<double> _poleMasses;
  std::vector< std::string> _gFactorNames;
  std::map<int, std::vector<double> > _gFactorMap;
  std::vector<std::shared_ptr<AbsPhaseSpace> > _phpVecs;
  std::shared_ptr<TMatrixRel> _tMatr;
  std::shared_ptr<KMatrixBase> _kMatr;
  std::vector< std::shared_ptr<KPole> > _kPoles;
  ParticleTable* _particleTable;
  
  TFile* _theTFile;
  std::vector<TH1F*> _AmpRealH1Vec;
  std::vector<TH1F*> _AmpImagH1Vec;
  std::vector<TH2F*> _ArgandH2Vec;
  std::vector<TH2F*> _PhaseH2Vec;
  std::vector<TH1F*> _ElasticityH1Vec;
  std::vector<TH1F*> _SqrT11H1Vec;
  std::vector<TH1F*> _SqrT1iH1Vec;
  std::vector<TH1F*> _T1iH1Vec;
  std::vector<TH1F*> _SqrS1iH1Vec;
  std::vector<TH1F*> _phpH1Vec;
  std::vector<TH1F*> _deltaiiVec; 
  std::vector<TH1F*> _delta1iVec;
  void init();
};


