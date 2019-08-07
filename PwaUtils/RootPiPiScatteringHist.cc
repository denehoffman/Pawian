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

#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <complex>
#include <boost/algorithm/string.hpp>

#include "PwaUtils/RootPiPiScatteringHist.hh"
#include "PwaUtils/EvtDataScatteringList.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/IdStringMapRegistry.hh"
#include "Utils/PawianConstants.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "ConfigParser/ParserBase.hh"
#include "FitParams/AbsPawianParameters.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TGraphErrors.h"
#include "TGraph.h"

#include "TLorentzVector.h"
#include "ErrLogger/ErrLogger.hh"
#include "TTree.h"

RootPiPiScatteringHist::RootPiPiScatteringHist(std::string additionalSuffix, bool withTruth) :
  AbsHist(additionalSuffix, withTruth)
 {
  std::ostringstream rootFileName;
  rootFileName << "./pawianHists" << GlobalEnv::instance()->outputFileNameSuffix() << _additionalSuffix.c_str() <<  ".root";
  _theTFile=new TFile(rootFileName.str().c_str(),"recreate");

  _phaseDataGraphErr = new TGraphErrors();
  _phaseDataGraphErr->SetName("deltaData");
  _phaseDataGraphErr->SetTitle("#delta_{ij} (data)");
  _etaDataGraphErr = new TGraphErrors();
  _etaDataGraphErr->SetName("etaData");
  _etaDataGraphErr->SetTitle("#eta_{ij} (data)");

  _argandUnitsDataGraphErr = new TGraphErrors();
  _argandUnitsDataGraphErr->SetName("ArgandUnitsData");
  _argandUnitsDataGraphErr->SetTitle("ArgandUnits (data)");


  _phaseFitGraphErr = new TGraphErrors();
  _phaseFitGraphErr->SetName("deltaFit");
  _phaseFitGraphErr->SetTitle("#delta_{ij} (fit)");
  _etaFitGraphErr = new TGraphErrors();
  _etaFitGraphErr->SetName("etaFit");
  _etaFitGraphErr->SetTitle("#eta_{ij} (fit)");
  _argandUnitsFitGraphErr = new TGraphErrors();
  _argandUnitsFitGraphErr->SetName("ArgandUnitsFit");
  _argandUnitsFitGraphErr->SetTitle("ArgandUnits (fit)");


  _dataFourvecs = new TTree("_dataFourvecs", "_dataFourvecs");
  _fittedFourvecs = new TTree("_fittedFourvecs", "_fittedFourvecs");

  _dataFourvecs->Branch("mass", &_massVal, "mass");
  _dataFourvecs->Branch("dataVal", &_dataVal, "dataVal");
  _dataFourvecs->Branch("dataErrVal", &_dataErrVal, "dataErrVal");

  _fittedFourvecs->Branch("mass", &_massVal, "mass");
  _fittedFourvecs->Branch("fitVal", &_fitVal, "fitVal");
  //  _fittedFourvecs->Branch("dataErrVal", &_dataErrVal, "dataErrVal");
 }

RootPiPiScatteringHist::~RootPiPiScatteringHist(){
  _phaseDataGraphErr->Write();
  _etaDataGraphErr->Write();
  _argandUnitsDataGraphErr->Write();

  _phaseFitGraphErr->Write();
  _etaFitGraphErr->Write();
  _argandUnitsFitGraphErr->Write();

  _theTFile->Write();
  _theTFile->Close();
}

void RootPiPiScatteringHist::fillFromLhData(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> fitParams){

  if(0==theLh){
    Alert <<"AbsLh* is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }
  double weight=1.; //dummy

  theLh->updateFitParams(fitParams);

  double weightFit = fitParams->Value(theLh->getChannelScaleParam());

  const std::vector<EvtData*> dataList=theLh->getDataVec();
  double integralDataWWeight=0.;

  std::vector<EvtData*>::const_iterator it=dataList.begin();
  int dataPoint=1;
  while(it!=dataList.end())
    {
      //      integralDataWWeight+=weight;
      fillEvt((*it), weight, "data", dataPoint);
      fillEvt((*it), weightFit, "fit", dataPoint);
      InfoMsg << "data No " << (*it)->evtNo << " filled!!!" << endmsg;
      ++dataPoint;
      ++it;
    }

}

void RootPiPiScatteringHist::fillEvt(EvtData* theData, double weight, std::string evtType, int pointNr){

  if(evtType=="mc") return; // no MC data available

  TTree* theTree=0;
  TGraphErrors* _dataGraph=0;

  if(evtType=="data"){
    theTree=_dataFourvecs;
    if(GlobalEnv::instance()->parser()->productionFormalism()=="Phase" || GlobalEnv::instance()->parser()->productionFormalism()=="PhaseDiff"){
    _dataGraph = _phaseDataGraphErr;
    }
    else if (GlobalEnv::instance()->parser()->productionFormalism()=="Elasticity"){
    _dataGraph =_etaDataGraphErr;
    }
    else if (GlobalEnv::instance()->parser()->productionFormalism()=="ArgandUnits"){
    _dataGraph =_argandUnitsDataGraphErr;
    }
  }

  else if(evtType=="fit"){
    theTree=_fittedFourvecs;
    if(GlobalEnv::instance()->parser()->productionFormalism()=="Phase" || GlobalEnv::instance()->parser()->productionFormalism()=="PhaseDiff"){
      _dataGraph = _phaseFitGraphErr;
    }
    else if (GlobalEnv::instance()->parser()->productionFormalism()=="Elasticity"){
    _dataGraph =_etaFitGraphErr;
    }
    else if (GlobalEnv::instance()->parser()->productionFormalism()=="ArgandUnits"){
      _dataGraph = _argandUnitsFitGraphErr;
    }
  }

  if(evtType=="data" || evtType=="fit"){
    _massVal=theData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::M_PIPISCAT_NAME)) ;
    if(evtType=="data"){
      _dataVal=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::DATA_PIPISCAT_NAME));
      _dataErrVal=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::DATAERR_PIPISCAT_NAME));
      _dataGraph->SetPoint(pointNr, _massVal, _dataVal);
      _dataGraph->SetPointError(pointNr, 0., _dataErrVal);
    }
    else{
      _fitVal=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME));
      _fitVal*=weight;
      _dataGraph->SetPoint(pointNr, _massVal, _fitVal);
      _dataGraph->SetPointError(pointNr, 0., 0.);
    }
    theTree->Fill();
  }
}



