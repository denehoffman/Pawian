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
//#include "PwaUtils/EvtDataBaseList.hh"
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
  _TijSqrDataGraphErr = new TGraphErrors();
  _TijSqrDataGraphErr->SetName("sqrTijData");
  _TijSqrDataGraphErr->SetTitle("|T_{ij}|^{2} (data)");
  _ArgandDataGraph = new TGraph();
  _ArgandDataGraph->SetName("ArgandPlotData");
  _ArgandDataGraph->SetTitle("Argand plot (data)");

  _phaseFitGraphErr = new TGraphErrors();
  _phaseFitGraphErr->SetName("deltaFit");
  _phaseFitGraphErr->SetTitle("#delta_{ij} (fit)");
  _etaFitGraphErr = new TGraphErrors();
  _etaFitGraphErr->SetName("etaFit");
  _etaFitGraphErr->SetTitle("#eta_{ij} (fit)");
  _TijSqrFitGraphErr = new TGraphErrors();
  _TijSqrFitGraphErr->SetName("sqrTijFit");
  _TijSqrFitGraphErr->SetTitle("|T_{ij}|^{2} (fit)");
  _ArgandFitGraph = new TGraph();
  _ArgandFitGraph->SetName("ArgandPlotFit");
  _ArgandFitGraph->SetTitle("Argand plot (fit)");  

  _dataFourvecs = new TTree("_dataFourvecs", "_dataFourvecs");
  _fittedFourvecs = new TTree("_fittedFourvecs", "_fittedFourvecs");

  _dataFourvecs->Branch("mass", &_massVal, "mass");
  _dataFourvecs->Branch("phi", &_phiVal, "phi");
  _dataFourvecs->Branch("phiErr", &_phiErrVal, "phiErr");
  _dataFourvecs->Branch("eta", &_etaVal, "eta");
  _dataFourvecs->Branch("etaErr", &_etaErrVal, "etaErr");

  _fittedFourvecs->Branch("mass", &_massVal, "mass");
  _fittedFourvecs->Branch("phi", &_phiVal, "phi");
  _fittedFourvecs->Branch("phiErr", &_phiErrVal, "phiErr");
  _fittedFourvecs->Branch("eta", &_etaVal, "eta");
  _fittedFourvecs->Branch("etaErr", &_etaErrVal, "etaErr");
 }

RootPiPiScatteringHist::~RootPiPiScatteringHist(){
  _phaseDataGraphErr->Write();
  _etaDataGraphErr->Write();
  _ArgandDataGraph->Write();
  _TijSqrDataGraphErr->Write();

  _phaseFitGraphErr->Write();
  _etaFitGraphErr->Write();
  _ArgandFitGraph->Write();
  _TijSqrFitGraphErr->Write();

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

  const std::vector<EvtData*> dataList=theLh->getDataVec();
  double integralDataWWeight=0.;

  std::vector<EvtData*>::const_iterator it=dataList.begin();
  int dataPoint=1;
  while(it!=dataList.end())
    {
      //      integralDataWWeight+=weight;
      fillEvt((*it), weight, "data", dataPoint);
      InfoMsg << "data No " << (*it)->evtNo << " filled!!!" << endmsg;
      ++dataPoint;
      ++it;
    }

}

void RootPiPiScatteringHist::fillEvt(EvtData* theData, double weight, std::string evtType, int pointNr){

  if(evtType=="mc") return; // no MC data available

  TTree* theTree=0;
  TGraphErrors* _phaseGraph=0;
  TGraphErrors* _etaGraph=0;
  TGraphErrors* _TijSqrGraphErr=0;
  TGraph* _ArgandGraph=0;

  if(evtType=="data"){
    theTree=_dataFourvecs;
    _phaseGraph=_phaseDataGraphErr;
    _etaGraph=_etaDataGraphErr;
    _TijSqrGraphErr=_TijSqrDataGraphErr;
    _ArgandGraph=_ArgandDataGraph;
  }
  else if(evtType=="fit"){
    theTree=_fittedFourvecs;
    _phaseGraph=_phaseFitGraphErr;
    _etaGraph=_etaFitGraphErr;
    _TijSqrGraphErr=_TijSqrFitGraphErr;
    _ArgandGraph=_ArgandFitGraph;
  }

  if(evtType=="data" || evtType=="fit"){
    int evtNo=pointNr;
    _massVal=theData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::M_PIPISCAT_NAME)) ;
    _phiVal=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::PHI_PIPISCAT_NAME));
    _phiErrVal=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::PHIERR_PIPISCAT_NAME));
    _etaVal=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::ETA_PIPISCAT_NAME));
    _etaErrVal=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::ETAERR_PIPISCAT_NAME));
    _phaseGraph->SetPoint(evtNo, _massVal, _phiVal);
    _phaseGraph->SetPointError(evtNo, 0., _phiErrVal);

    _etaGraph->SetPoint(evtNo, _massVal, _etaVal);
    _etaGraph->SetPointError(evtNo, 0., _etaErrVal);

    double twoDeltaRad=2.*_phiVal*PawianConstants::degToRad-PawianConstants::pi/2.;
    double mag=_etaVal/2.;
    complex<double> argendPoint= complex<double>(mag*cos(twoDeltaRad), mag*sin(twoDeltaRad)+0.5);
    _ArgandGraph->SetPoint(evtNo, argendPoint.real(), argendPoint.imag());

    complex<double> etaVec2EiDelta11(_etaVal*cos(2.*_phiVal*PawianConstants::degToRad), _etaVal*sin(2.*_phiVal*PawianConstants::degToRad));
    complex<double> currentT11 = (etaVec2EiDelta11-complex<double>(1.,0.))/(2.*PawianConstants::i);
    _TijSqrGraphErr->SetPoint(evtNo, _massVal, norm(currentT11));
    theTree->Fill();
  }
}



