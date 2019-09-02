//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// RootPiPiScatteringHist class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <sstream>
#include <memory>

#include "TROOT.h"

#include "PwaUtils/AbsHist.hh"
#include "Utils/PawianCollectionUtils.hh"

class AbsLh;
class EvtData;
class TLorentzVector; 

class TFile;
class TGraphErrors;
class TGraph;
class TTree;
class Particle;
class AbsPawianParameters;

class RootPiPiScatteringHist : public AbsHist{

public:
  RootPiPiScatteringHist(std::string additionalSuffix="", bool withTruth=false);
  virtual ~RootPiPiScatteringHist();
  virtual void fillEvt(EvtData* theData, double weight, std::string evtType, int pointNr);
  virtual void fillFromLhData(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> fitParams);
  virtual void scaleFitHists(double scaleFactor) {;}
  
protected:
  TFile* _theTFile;
  TTree* _dataFourvecs;
  TTree* _fittedFourvecs;

  TGraphErrors* _dataGraphErr;
  TGraphErrors* _fitGraphErr;

  //  virtual void initRootStuff()=0;
  
private:
  float _massVal;
  float _dataVal;
  float _dataErrVal;
  float _fitVal;
  float _fitErrVal;
};

