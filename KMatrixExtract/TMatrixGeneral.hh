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

#include "KMatrixExtract/TMatrixGeneralBase.hh"

class TFile;
class TH1F;
class TH2F;

class TMatrixGeneral : public TMatrixGeneralBase {

public:

  // create/copy/destroy:

  ///Constructor
  TMatrixGeneral(); 
  TMatrixGeneral(pipiScatteringParser* theParser);


  /** Destructor */
  virtual ~TMatrixGeneral();

  virtual void process();
  virtual void initHistos();

  // Getters:
 
protected:
  unsigned int _noOfSteps;
  double _stepSize; 
  std::vector<double> _energyPlaneBorders;
  unsigned int _numXStepsForSheetScan;
  unsigned int _numYStepsForSheetScan;
  TFile* _theTFile;
  bool _histosInit;

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
};


