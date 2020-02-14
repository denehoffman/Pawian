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
#include "pipiScatteringUtils/FVectorIntensityGeneral.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh" 
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/TMatrixBase.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/KMatrixBase.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/KPoleBarrier.hh"
#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/KMatrixRelBg.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceFactory.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaUtils/FVectorIntensityDynamics.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/AbsDecay.hh"
#include "pipiScatteringUtils/PiPiScatteringChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ConfigParser/KMatrixParser.hh"
#include "ConfigParser/pipiScatteringParser.hh"

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
//#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

FVectorIntensityGeneral::FVectorIntensityGeneral(pipiScatteringParser* theParser) :
  _pipiScatteringParser(theParser)
  ,_pipiScatteringChannelEnv(new PiPiScatteringChannelEnv(theParser))
  ,_projectionParticleNames("")
  ,_motherParticleName("")
  ,_pVecName("")
  ,_decProjectionIndex(0)
  ,_noOfSteps(1000)
  ,_massMin(1.)
  ,_massMax(2.)
  ,_stepSize(0.001)
  ,_pathToFitParams("") 
{ 
  init();
  std::string rootFileName="./FVectorIntensityGeneral.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");

  std::string magKey="Magnitude "+ _pVecName + " to " +_projectionParticleNames;
  std::string phaseKey="Phase "+ _pVecName + " to " +_projectionParticleNames;
  std::string intensityKey="Intensity "+ _pVecName + " to " +_projectionParticleNames;

  _MagH1 = new TH1F(magKey.c_str(), magKey.c_str(), 
                                  _noOfSteps-1, _massMin, _massMax);

  _PhaseH1 = new TH1F(phaseKey.c_str(), phaseKey.c_str(), 
                                  _noOfSteps-1, _massMin, _massMax);

  _IntensityH1 = new TH1F(intensityKey.c_str(), intensityKey.c_str(), 
                                  _noOfSteps-1, _massMin, _massMax);
}

FVectorIntensityGeneral::~FVectorIntensityGeneral() {
  _theTFile->Write();
  _theTFile->Close();
}


void FVectorIntensityGeneral::init() {
  GlobalEnv::instance()->setup();
  _particleTable=GlobalEnv::instance()->particleTable();
  GlobalEnv::instance()->AddEnv(_pipiScatteringChannelEnv, AbsChannelEnv::CHANNEL_PIPISCATTERING);
  _pipiScatteringChannelEnv->setupChannel(0);

  std::string pathToKMatrixParser=_pipiScatteringChannelEnv->pathToKMatrixParser();
  InfoMsg << "pathToKMatrixParser: " << pathToKMatrixParser << endmsg;

  _pathToFitParams = _pipiScatteringParser->fitParamFile();
  InfoMsg << "path th fit parameters: " << _pathToFitParams << endmsg;

  std::vector<Particle*> fsParticles = _pipiScatteringChannelEnv->finalStateParticles();
  std::vector<Particle*>::iterator it;
  for(it=fsParticles.begin(); it!=fsParticles.end(); ++it){
    InfoMsg << "fsParticle: " << (*it)->name() << endmsg;
  }

  std::shared_ptr<AbsDecayList> absDecList=_pipiScatteringChannelEnv->absDecayList();
  std::vector<std::shared_ptr<AbsDecay> > theDecs=absDecList->getList();
  if(theDecs.size()!=1){
    Alert << "the decay list contains " << theDecs.size() << " decays" << endmsg;
    Alert << "exactly one decay is required!!! " << endmsg;
    exit(1);    
  }

  std::shared_ptr<AbsDecay> theDec=theDecs.at(0);
  Particle* theMotherParticle=theDec->motherPart();
  _motherParticleName=theMotherParticle->name();
  InfoMsg << "theMotherParticle: " << _motherParticleName << endmsg;

  _projectionParticleNames= theDec->projectionParticleNames();   
  InfoMsg << "projectionParticleNames: " << _projectionParticleNames << endmsg;

  std::string baseNameFVector=_pipiScatteringParser->baseNameFVector();
  InfoMsg << "baseNameFVector: " << baseNameFVector << endmsg;

 
  ChannelID channelID(0);
  std::string decayDynamics=_pipiScatteringParser->decayDynamics().at(0);
  InfoMsg << "decayDynamics: " << decayDynamics << endmsg; 

  _pVecName=baseNameFVector+"b"+_motherParticleName;
  
   std::string dummyName="dummy"; 
  _fVectorIntensityDynamics = 
    std::shared_ptr<FVectorIntensityDynamics>(new FVectorIntensityDynamics(dummyName, fsParticles, theMotherParticle, pathToKMatrixParser, baseNameFVector, channelID, _projectionParticleNames));
  _orbitalL= _fVectorIntensityDynamics->orbitalL();
  
  _phpVecs = _fVectorIntensityDynamics->getKMatix()->phaseSpaceVec();
  _decProjectionIndex = _fVectorIntensityDynamics->decProjectionIndex();
  InfoMsg << "_decProjectionIndex: " << _decProjectionIndex << endmsg;

  _phpVecCurrent = _phpVecs.at(_decProjectionIndex);
  _massMin = _phpVecCurrent->thresholdMass()+0.001;
  _stepSize=(_massMax-_massMin)/_noOfSteps;

  std::shared_ptr<AbsPawianParameters> params=ParamFactory::instance()->getParametersPointer("Pawian");
  _fVectorIntensityDynamics->fillDefaultParams(params);

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
  if(_pathToFitParams != "") _fVectorIntensityDynamics->updateFitParams(_params);

  _fVector=_fVectorIntensityDynamics->getFVector();  
}

void FVectorIntensityGeneral::process(){
  for (double mass=_massMin+_stepSize/0.5; mass<_massMax; mass+=_stepSize){
    complex<double> currentResult=_fVector->evalProjMatrix(mass, _decProjectionIndex, _orbitalL);
    _MagH1->Fill(mass, std::abs(currentResult));
    _PhaseH1->Fill(mass, std::arg(currentResult)*PawianConstants::radToDeg);
    _IntensityH1->Fill(mass, norm( currentResult*sqrt(_phpVecCurrent->factor(mass).real())) );
  }
}

