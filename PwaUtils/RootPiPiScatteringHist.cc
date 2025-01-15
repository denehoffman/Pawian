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

  std::string graphNameData("DeltaData");
  std::string graphNameFit("DeltaFit");

  std::string graphTitleData="#delta^{o} (data)";
  std::string graphTitleFit("#delta^{o} (fit)");

  if (GlobalEnv::instance()->parser()->productionFormalism()=="Phase"){
    graphNameData="PhaseData";
    graphNameFit="PhaseFit";
    graphTitleData="#delta (data)";
    graphTitleFit="#delta (fit)";
  }
  else if (GlobalEnv::instance()->parser()->productionFormalism()=="PhasePhi"){
    graphNameData="PhaseData";
    graphNameFit="PhaseFit";
    graphTitleData="#Phi (data)";
    graphTitleFit="#Phi (fit)";
  }
  else if (GlobalEnv::instance()->parser()->productionFormalism()=="Elasticity"){
    graphNameData="ElasticityData";
    graphNameFit="ElasticityFit";
    graphTitleData="#eta (data)";
    graphTitleFit="#eta (fit)";
  }
  else if (GlobalEnv::instance()->parser()->productionFormalism()=="ArgandUnits"){
    graphNameData="ArgandUnitsData";
    graphNameFit="ArgandUnitsFit";
    graphTitleData="#rho_{i} #rho_{j} |T_{ij}|^{2} (data)";                                
    graphTitleFit="#rho_{i} #rho_{j} |T_{ij}|^{2} (fit)"; 
  }
  else if (GlobalEnv::instance()->parser()->productionFormalism()=="Treal"){
    graphNameData="TrealData";
    graphNameFit="TrealFit";
    graphTitleData="T_{ij} real (data)";
    graphTitleFit="T_{ij} real (fit)";
  }
  else if (GlobalEnv::instance()->parser()->productionFormalism()=="Timag"){
    graphNameData="TimagData";
    graphNameFit="TimagFit";
    graphTitleData="T_{ij} imag (data)";
    graphTitleFit="T_{ij} imag (fit)";
  }
  else if (GlobalEnv::instance()->parser()->productionFormalism()=="Tabs"){
    graphNameData="TabsData";
    graphNameFit="TabsFit";
    graphTitleData="|T| (data)";
    graphTitleFit="|T| (fit)";
  }
  else if (GlobalEnv::instance()->parser()->productionFormalism()=="PhaseDiff"){
    graphNameData="PhaseDiffData";
    graphNameFit="PhaseDiffFit";
    graphTitleData="#phi_{i} - #phi_{j}  (data)";
    graphTitleFit="#phi_{i} - #phi_{j} (fit)";
  }
  else if (GlobalEnv::instance()->parser()->productionFormalism()=="PVecIntensity"){
    graphNameData="PVecIntensityData";
    graphNameFit="PVecIntensityFit";
    graphTitleData="intensity (data)";
    graphTitleFit="intensity (fit)";
  }
  else{
    Alert <<"production formalism " 
	  << GlobalEnv::instance()->parser()->productionFormalism() 
	  << " is not supported here!!!" << endmsg;
    exit(1);
  }


  _dataGraphErr = new TGraphErrors();
  _dataGraphErr->SetName(graphNameData.c_str());
  _dataGraphErr->SetTitle(graphTitleData.c_str());

  _fitGraphErr = new TGraphErrors();
  _fitGraphErr->SetName(graphNameFit.c_str());
  _fitGraphErr->SetTitle(graphTitleFit.c_str());

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
  _dataGraphErr->Write();
  _fitGraphErr->Write();

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

  std::vector<EvtData*>::const_iterator it=dataList.begin();
  int dataPoint=1;
  while(it!=dataList.end())
    {
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

  if(evtType=="data" || evtType=="fit"){
    _massVal=theData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::M_PIPISCAT_NAME)) ;
    if(evtType=="data"){
      theTree=_dataFourvecs;
      _dataVal=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::DATA_PIPISCAT_NAME));
      _dataErrVal=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::DATAERR_PIPISCAT_NAME));
	_dataGraphErr->SetPoint(pointNr, _massVal, _dataVal);
	_dataGraphErr->SetPointError(pointNr, 0., _dataErrVal);
    }
    else{  //evtType=="fit"
      theTree=_fittedFourvecs;
      _fitVal=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME));
      _fitErrVal=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FITERR_PIPISCAT_NAME));
      _fitVal*=weight;
      _fitGraphErr->SetPoint(pointNr, _massVal, _fitVal);
      _fitGraphErr->SetPointError(pointNr, 0., _fitErrVal);
    }
    theTree->Fill();
  }
}



