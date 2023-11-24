//************************************************************************//
//									  //
//  Copyright 2019 Bertram Kopf (bertram@ep1.rub.de)			  //
//          	   - Ruhr-Universit??t Bochum 				  //
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

// TMatrixCompareDynamics class definition file. -*- C++ -*-
// Copyright 2019 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "Utils/PawianIOUtils.hh"
#include "PwaUtils/TMatrixCompareDynamics.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataScatteringList.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "ConfigParser/KMatrixParser.hh"
#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/KMatrixRelBg.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/KPoleBarrier.hh"
#include "PwaDynamics/PPole.hh"
#include "PwaDynamics/PPoleBarrier.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceFactory.hh"
#include "PwaDynamics/KMatrixFunctions.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "Utils/IdStringMapRegistry.hh"
#include "Utils/PawianConstants.hh"

TMatrixCompareDynamics::TMatrixCompareDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser,  std::string& pathToKMatrCompareConfigFile, std::string dataType, std::string projectionParticleNames) :
  TMatrixDynamics(name, fsParticles, mother, pathToConfigParser, dataType, projectionParticleNames)
  ,_projectionCompareIndex(0)
  ,_currentOffset(0.)
{
  std::string completePathToComparCfgConfig=PawianIOUtils::getFileName(GlobalEnv::instance()->KMatrixStorePath(),pathToKMatrCompareConfigFile); 
  
  _kMatrixParserCompare= std::shared_ptr<KMatrixParser>(new KMatrixParser(completePathToComparCfgConfig));
  _tMatrDynCompare=std::shared_ptr<TMatrixDynamics>(new TMatrixDynamics(_kMatrixParserCompare));
  _tMatrCompare=_tMatrDynCompare->getTMatix();  
  _orbitalLCompare=_kMatrixParserCompare->orbitalMom();
  
  std::istringstream projParticles(_projectionParticleNames);
  std::string firstProjParticleName;
  std::string secondProjParticleName;
  projParticles >> firstProjParticleName >> secondProjParticleName;
  std::string projKey=firstProjParticleName+secondProjParticleName;    

  bool found=false;
  std::vector< std::string>   gFactorNamesCompare=_tMatrDynCompare->gFactorNames();
    for(unsigned int idx=0; idx<gFactorNamesCompare.size();++idx){
    if(projKey==gFactorNamesCompare.at(idx)){
      _projectionCompareIndex=idx;
      found=true;
    }
  }
  if (!found){
    Alert << "projection index for key " << projKey << " not found in TMatrix compare" << endmsg;
    exit(0);
  }
  InfoMsg << _tMatrDynCompare->KMatixName() << ": projection index: " << _projectionCompareIndex <<endmsg;
}

TMatrixCompareDynamics::~TMatrixCompareDynamics()
{
}

complex<double> TMatrixCompareDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=_tMatr->kMatrix()->phaseSpaceVec();
  double currentMass=theData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::M_PIPISCAT_NAME));
  _tMatr->evalMatrix(currentMass, _orbitalL);
  _tMatrCompare->evalMatrix(currentMass, _orbitalLCompare);
  evalPhaseCompare(theData, currentMass);
  return (*_tMatr)(_prodProjectionIndex,_decProjectionIndex);
}

void TMatrixCompareDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  TMatrixDynamics::fillDefaultParams(fitPar);
  _tMatrDynCompare->fillDefaultParams(fitPar);
  std::string offsetName=name()+_tMatrDynCompare->name()+"PhaseDiffOffset";
  fitPar->Add(offsetName, 0., 0.1);
  fitPar->SetLimits( offsetName, -360., 360.);
}

void TMatrixCompareDynamics::fillParamNameList(){
  _paramNameList.clear();
  TMatrixDynamics::fillParamNameList();
  _tMatrDynCompare->fillParamNameList();
  std::string offsetName=name()+_tMatrDynCompare->name()+"PhaseDiffOffset";
  _paramNameList.push_back(offsetName);
}

bool TMatrixCompareDynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){
  return true;
}

void TMatrixCompareDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  TMatrixDynamics::updateFitParams(fitPar);
  _tMatrDynCompare->updateFitParams(fitPar);
  std::string offsetName=name()+_tMatrDynCompare->name()+"PhaseDiffOffset";
  _currentOffset=fitPar->Value(offsetName);
}

void TMatrixCompareDynamics::fillMasses(EvtData* theData){
  TMatrixDynamics::fillMasses(theData);
  _tMatrDynCompare->fillMasses(theData);
}

void TMatrixCompareDynamics::evalPhaseCompare(EvtData* theData, double currentMass){
  double currentPhase=std::arg((*_tMatr)(_decProjectionIndex,_decProjectionIndex));
  double currentPhaseCompare=std::arg((*_tMatrCompare)(_decProjectionIndex,_decProjectionIndex));
  double currentPhaseDiff=(currentPhase-currentPhaseCompare)*PawianConstants::radToDeg;
  while(currentPhaseDiff > 180. ) currentPhaseDiff-=180.;
  while(currentPhaseDiff < -180.) currentPhaseDiff+=180.;
  theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME))=currentPhaseDiff+_currentOffset;
 }


