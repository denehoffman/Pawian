//************************************************************************//
//									  //
//  Copyright 2023 Bertram Kopf (bertram@ep1.rub.de)			  //
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
#include "KMatrixExtract/TMatrixGeneralBase.hh"
#include "KMatrixExtract/RiemannSheetAnalyzer.hh"
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
#include "PwaDynamics/KMatrixFunctions.hh"
#include "ConfigParser/KMatrixParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/PdtParser.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "FitParams/AbsPawianParamStreamer.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "Utils/PawianConstants.hh"
#include "pipiScatteringUtils/PiPiScatteringChannelEnv.hh"
#include "ConfigParser/pipiScatteringParser.hh"

#include "ConfigParser/ParserBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/TMatrixDynamics.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "FitParams/ParamFactory.hh"

#include "ErrLogger/ErrLogger.hh"

TMatrixGeneralBase::TMatrixGeneralBase(pipiScatteringParser* theParser) :
  _pipiScatteringParser(theParser)
  ,_pipiScatteringChannelEnv(new PiPiScatteringChannelEnv(theParser))
  ,_pathToFitParams("")
  ,_orbitalL(0)
  ,_pathToKMatrixParser("")
  ,_massMin(100000.)
  ,_massMax(0.) 
{
  init();
}

TMatrixGeneralBase::TMatrixGeneralBase(){
}

TMatrixGeneralBase::~TMatrixGeneralBase() {
}


void TMatrixGeneralBase::init() {
  GlobalEnv::instance()->setup();
  _particleTable=GlobalEnv::instance()->particleTable();
  GlobalEnv::instance()->AddEnv(_pipiScatteringChannelEnv, AbsChannelEnv::CHANNEL_PIPISCATTERING);
  _pipiScatteringChannelEnv->setupChannel(0);
  
  _pathToKMatrixParser=_pipiScatteringChannelEnv->pathToKMatrixParser();
  InfoMsg << "pathToKMatrixParser: " << _pathToKMatrixParser << endmsg;

  _pathToFitParams = _pipiScatteringParser->fitParamFile();
  InfoMsg << "path th fit parameters: " << _pathToFitParams << endmsg;

  _fsParticles = _pipiScatteringChannelEnv->finalStateParticles();
  std::vector<Particle*>::iterator it;
  for(it=_fsParticles.begin(); it!=_fsParticles.end(); ++it){
    InfoMsg << "fsParticle: " << (*it)->name() << endmsg;
  }

  InfoMsg << "_fsParticles.size(): " << _fsParticles.size() << endmsg;

  std::shared_ptr<AbsDecayList> absDecList=_pipiScatteringChannelEnv->absDecayList();
  std::vector<std::shared_ptr<AbsDecay> > theDecs=absDecList->getList();
  if(theDecs.size()!=1){
    Alert << "the decay list contains " << theDecs.size() << " decays" << endmsg;
    Alert << "exactly one decay is required!!! " << endmsg;
    exit(1);    
  }

  std::shared_ptr<AbsDecay> theDec=theDecs.at(0);
  _motherParticle=theDec->motherPart();
  _motherParticleName=_motherParticle->name();
  InfoMsg << "theMotherParticle: " << _motherParticleName << endmsg;

  _projectionParticleNames= theDec->projectionParticleNames();   
  InfoMsg << "projectionParticleNames: " << _projectionParticleNames << endmsg;



  _kMatrixParser= std::shared_ptr<KMatrixParser>(new KMatrixParser(_pathToKMatrixParser));

  _tMatrDyn = std::shared_ptr<TMatrixDynamics>(new TMatrixDynamics(_kMatrixParser));

  _kMatr = _tMatrDyn->getKMatix(); 
  _tMatr =  _tMatrDyn->getTMatix();  

  _phpVecs=_kMatr->phaseSpaceVec();
  _gFactorNames= _tMatrDyn->gFactorNames();

  _orbitalL= _tMatrDyn->orbitalL();
  _decProjectionIndex = _tMatrDyn->decProjectionIndex();

  std::vector<std::string> poleNameAndMassVecs=_kMatrixParser->poles();
  std::vector<std::string>::iterator itString;
  for (itString=poleNameAndMassVecs.begin(); itString!=poleNameAndMassVecs.end(); ++itString){
    std::istringstream poleIString(*itString);
    std::string currentPoleName;
    std::string currentPoleMassStr;
    poleIString >> currentPoleName >> currentPoleMassStr;
    
    std::istringstream currentPoleMassiStr(currentPoleMassStr);
    double currentValue;
    if(!(currentPoleMassiStr >> currentValue)){
      Alert << "cannot convert " << currentPoleMassStr << " to a double value" << endmsg;
      exit(0);
    }
    if (currentValue>_massMax) _massMax=currentValue;
    
  }

  const std::vector<std::string> gFacStringVec=_kMatrixParser->gFactors();
  DebugMsg << "gFacStringVec.size(): " << gFacStringVec.size() << endmsg;
  std::map<std::pair<std::string, std::string>, std::string> phpDescriptionVec =
    _kMatrixParser->phpDescriptionMap();
  std::cout << "phpDescriptionVec.size(): " << phpDescriptionVec.size() << std::endl;
  std::vector<std::string>::const_iterator itStrConst;
  for(itStrConst=gFacStringVec.begin(); itStrConst!=gFacStringVec.end(); ++itStrConst){
    std::istringstream particles(*itStrConst);
    std::string firstParticleName;
    std::string secondParticleName;
    particles >> firstParticleName >> secondParticleName;
    std::pair<std::string, std::string> currentParticlePair=make_pair(firstParticleName, secondParticleName);
    Particle* firstParticle = _particleTable->particle(firstParticleName);
    Particle* secondParticle = _particleTable->particle(secondParticleName);
    if(0==firstParticle || 0==secondParticle){
      Alert << "particle with name: " << firstParticleName <<" or " << secondParticleName 
	    << " doesn't exist in pdg-table" << endmsg;
      exit(0);
    }
    double currentMassSum=firstParticle->mass()+secondParticle->mass();
    if(currentMassSum<_massMin) _massMin=currentMassSum;
  }

}

void TMatrixGeneralBase::fillParams(){

  std::shared_ptr<AbsPawianParameters> params=ParamFactory::instance()->getParametersPointer("Pawian");
  _tMatrDyn->fillDefaultParams(params);
  _tMatrDyn->fillParamNameList();
 
    std::ifstream ifs(_pathToFitParams);
    if(!ifs.good()) 
      { //file doesn't exist; dum default params
	WarningMsg << "could not parse " << _pathToFitParams << endmsg;
	WarningMsg << "dump default parameter " << _pathToFitParams << endmsg;
	std::string defaultparamsname="defaultParams.dat";
	std::ofstream theStreamDefault ( defaultparamsname );
	params->print(theStreamDefault);
	theStreamDefault.close();
	exit(1);        
      }   


  AbsPawianParamStreamer thePawianStreamer(_pathToFitParams);
  _params = thePawianStreamer.paramList();

  InfoMsg << "The k-Matrix input parameter are: " << endmsg;
  _params->print(std::cout);
  if(_pathToFitParams != "") _tMatrDyn->updateFitParams(_params);

  _kMatrixParamNames = _tMatrDyn->paramNames();
  InfoMsg << "_kMatrixParamNames.size(): " << _kMatrixParamNames.size() << endmsg; 
}


void TMatrixGeneralBase::SetParamValue(const std::string & paramName, double paramVal){
  InfoMsg << "set param: " << paramName << " to " << paramVal << endmsg;
  _params->SetValue(paramName, paramVal);
}

double TMatrixGeneralBase::GetParamValue(const std::string & paramName){
  return _params->Value(paramName);
}

