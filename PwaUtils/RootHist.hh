//************************************************************************//
//									  //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// AbsHist class definition file. -*- C++ -*-
// Copyright 2016 Bertram Kopf

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
class TH2F;
class TH1F;
class TTree;
class Particle;
class AbsPawianParameters;

class RootHist : public AbsHist{

public:
  RootHist(std::string additionalSuffix="", bool withTruth=false);
  virtual ~RootHist();
  virtual void fillEvt(EvtData* theData, double weight, std::string evtType, int pointNr);
  virtual void fillFromLhData(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> fitParams);
  virtual void scaleFitHists(double scaleFactor);
  
protected:
  
  void fillMassHists(EvtData* theData, double weight, std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >& toFill);
  void fillAngleHists(EvtData* theData, double weight, std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess >& toFill, std::string frame="heli");
  void fillAngleHists2D(EvtData* theData, double weight, std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess >& toFill);

  TFile* _theTFile;
  TTree* _dataFourvecs;
  TTree* _fittedFourvecs;
  TTree* _truthFourvecs;

  bool _fillTruths;
  std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess > _massDataHistMap;
  std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess > _massMcHistMap;
  std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess > _massFitHistMap;
  std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess > _massTruthHistMap;
  std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess > _massTruthFitHistMap;
  
  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleDataHistMap;
  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleMcHistMap;
  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleFitHistMap;
  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleTruthHistMap;
  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleTruthFitHistMap;
  
  std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess > _angleDataHistMap2D;
  std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess > _angleMcHistMap2D;
  std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess > _angleFitHistMap2D;
  // std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess > _angleTruthHistMap2D;
  //   std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess > _angleTruthFitHistMap2D;
  
  
  //angles in the Gottfried Jackson system
  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleGJDataHistMap;
  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleGJMcHistMap;
  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleGJFitHistMap;
  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleGJTruthHistMap;
  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleGJTruthFitHistMap;
  
  //  virtual void initRootStuff()=0;
  
private:
  float _weightToWrite;
  std::map<std::string, std::shared_ptr<TLorentzVector> > _fourVecMap;
  
  virtual void initAngleHists(std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess >& theMap, std::shared_ptr<angleHistData> theHistData, std::string dataType, std::string systemType="Heli");
};

