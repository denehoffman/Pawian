//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "KMatrixExtract/FVectorIntensityGeneral.hh"
#include "KMatrixExtract/RiemannSheetFVectorAnalyzer.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh" 
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/TMatrixBase.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/KMatrixBase.hh"
#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/KMatrixRelBg.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaUtils/FVectorIntensityDynamics.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/AbsDecay.hh"
#include "pipiScatteringUtils/PiPiScatteringChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ConfigParser/KMatrixParser.hh"
#include "ConfigParser/pipiScatteringParser.hh"
#include "pipiScatteringUtils/PiPiScatteringChannelEnv.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Particle/PdtParser.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "FitParams/AbsPawianParamStreamer.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "Utils/PawianConstants.hh"

#include "ConfigParser/ParserBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/TMatrixDynamics.hh"
#include "FitParams/ParamFactory.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"


#include "ErrLogger/ErrLogger.hh"

FVectorIntensityGeneral::FVectorIntensityGeneral(pipiScatteringParser* theParser) :
  TMatrixGeneral(theParser)
  ,_pVecName("")
{ 
  init();
}

void FVectorIntensityGeneral::initHistos(){
  _histosInit=true;
  std::string rootFileName="./FVectorIntensityGeneral.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");

  std::string magKey="Magnitude "+ _pVecName + " to " + _gFactorNames.at(_decProjectionIndex);
  std::string phaseKey="Phase "+ _pVecName + " to " + _gFactorNames.at(_decProjectionIndex);
  std::string intensityKey="Intensity "+ _pVecName + " to " + _gFactorNames.at(_decProjectionIndex);

  std::string currentKey;
  for(unsigned int i=0; i < _gFactorNames.size(); ++i){
    currentKey="Magnitude"+ _pVecName + "to" + _gFactorNames.at(i);
    TH1F* currentMagH1 = new TH1F(currentKey.c_str(), currentKey.c_str(), _noOfSteps-1, _massMin, _massMax);
    _MagsH1.push_back(currentMagH1);

    currentKey="Phase"+ _pVecName + "to" + _gFactorNames.at(i);
    TH1F* currentPhaseH1 = new TH1F(currentKey.c_str(), currentKey.c_str(), 
                                  _noOfSteps-1, _massMin, _massMax);
    _PhasesH1.push_back(currentPhaseH1);

    currentKey="Intensity"+ _pVecName + "to" + _gFactorNames.at(i);
    TH1F* currentIntensityH1 = new TH1F(currentKey.c_str(), currentKey.c_str(),_noOfSteps-1, _massMin, _massMax);
    _IntensitiesH1.push_back(currentIntensityH1);

    currentKey="Argand"+ _pVecName + "to" + _gFactorNames.at(i);
    TGraph* currentArgandTGraph = new TGraph();
    currentArgandTGraph->SetName(currentKey.c_str());
    currentArgandTGraph->SetTitle(currentKey.c_str());

    _ArgandPlotsTGraph.push_back(currentArgandTGraph);
    
  }
}

FVectorIntensityGeneral::~FVectorIntensityGeneral() {
  if(_histosInit){
    for(unsigned int i=0; i < _gFactorNames.size(); ++i){  
      _ArgandPlotsTGraph.at(i)->Write();
    }
  }
}


void FVectorIntensityGeneral::init() {

  std::string baseNameFVector=_pipiScatteringParser->baseNameFVector();
  InfoMsg << "baseNameFVector: " << baseNameFVector << endmsg;
  std::vector<std::string> prodSuffixes= _pipiScatteringParser->addKmatrixProdSuffix();
 
  ChannelID channelID(0);
  //  _pVecName=baseNameFVector+"b"+_motherParticleName;

   std::string dummyName="dummy";

   InfoMsg << "_fsParticles.size(): " << _fsParticles.size() << endmsg;
   InfoMsg << "_pathToKMatrixParser: " << _pathToKMatrixParser << endmsg;
   InfoMsg << "baseNameFVector: " << baseNameFVector << endmsg;
   InfoMsg << "channelID: " << channelID << endmsg;
   InfoMsg << "_projectionParticleNames: " << _projectionParticleNames << endmsg;
   InfoMsg << "_motherParticle->name(): " << _motherParticle->name() << endmsg;
 
  _fVectorIntensityDynamics = 
    std::shared_ptr<FVectorIntensityDynamics>(new FVectorIntensityDynamics(dummyName, _fsParticles, _motherParticle, _pathToKMatrixParser, baseNameFVector, channelID, _projectionParticleNames));
  
  _decProjectionIndex = _fVectorIntensityDynamics->decProjectionIndex();
  InfoMsg << "_decProjectionIndex: " << _decProjectionIndex << endmsg;

  _phpVecCurrent = _phpVecs.at(_decProjectionIndex);

  fillParams();  
  _fVector=_fVectorIntensityDynamics->getFVector();
}

void FVectorIntensityGeneral::fillParams(){
  std::shared_ptr<AbsPawianParameters> params=ParamFactory::instance()->getParametersPointer("Pawian");
  //_fVectorIntensityDynamics->fillDefaultParams(params);

  std::ifstream ifs(_pathToFitParams);
  if(!ifs.good()) 
    { //file doesn't exist; dum default params
      WarningMsg << "could not parse " << _pathToFitParams << endmsg;
      WarningMsg << "dump defaut parameter " << _pathToFitParams << endmsg;
      std::string defaultparamsname="defaultParams.dat";
      std::ofstream theStreamDefault ( defaultparamsname );
      params->print(theStreamDefault);
      theStreamDefault.close();
      exit(1);        
    }

  AbsPawianParamStreamer thePawianStreamer(_pathToFitParams);
  _params = thePawianStreamer.paramList();

  InfoMsg << "The F-Vector input parameter are: " << endmsg;
  _params->print(std::cout);

  _fVectorIntensityDynamics->fillParamNameList();
  //  if(_pathToFitParams != "") _fVectorIntensityDynamics->updateFitParams(_params);
  std::vector<std::string> fVecParamNames=_fVectorIntensityDynamics->paramNames();
  for (unsigned int i=0; i<fVecParamNames.size(); ++i){
    InfoMsg << "F-Vector ParamNames: " << fVecParamNames.at(i) << endmsg;
  }
  _fVectorIntensityDynamics->fillDefaultParams(params);
  //  _fVectorIntensityDynamics->updateFitParams(params);
  InfoMsg << "The F-Vector input paras are: " << endmsg;
  params->print(std::cout);
}

void FVectorIntensityGeneral::process(){
  _fVectorIntensityDynamics->updateFitParams(_params);
  std::vector<int> pointNr;
  pointNr.resize(_gFactorNames.size());

  for (double mass=_massMin+_stepSize/0.5; mass<_massMax; mass+=_stepSize){
    _fVector->evalMatrix(mass, _orbitalL);
    for(unsigned int i=0; i < _gFactorNames.size(); ++i){
	complex<double> currentResult = (*_fVector)(i,0);
	_MagsH1.at(i)->Fill(mass, std::abs(currentResult));
	_PhasesH1.at(i)->Fill(mass, std::arg(currentResult)*PawianConstants::radToDeg);
	if(mass>_phpVecs.at(i)->thresholdMass()){
	  _IntensitiesH1.at(i)->Fill(mass, norm( currentResult*sqrt(_phpVecs.at(i)->factor(mass).real())));
	  _ArgandPlotsTGraph.at(i)->SetPoint(pointNr[i], currentResult.real(), currentResult.imag());
	  ++pointNr[i];
	} 
    }
  }
  RiemannSheetFVectorAnalyzer(_kMatrixParser->noOfChannels(), _fVector,
                       std::complex<double>(_energyPlaneBorders[0], _energyPlaneBorders[1]),
                       std::complex<double>(_energyPlaneBorders[2], _energyPlaneBorders[3]),
                       _numXStepsForSheetScan,
		       _numYStepsForSheetScan,
                       _decProjectionIndex);
}

