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

// FVectorCompareDynamics class definition file. -*- C++ -*-
// Copyright 2019 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/FVectorCompareDynamics.hh"
#include "PwaUtils/KMatrixDynamics.hh"
#include "PwaDynamics/FVector.hh"
#include "ConfigParser/KMatrixParser.hh"

#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataScatteringList.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
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

FVectorCompareDynamics::FVectorCompareDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother1, Particle* mother2, std::string& pathToConfigParser, std::string& pathToConfigParserComp, std::string dataType, ChannelID channelId, std::string projectionParticleNames) :
KMatrixDynamics(name, fsParticles, mother1, pathToConfigParser, channelId, projectionParticleNames)
  ,_projectionCompareIndex(0)
  , _kMatCompareName(name +"Compare")
  , _nameOfFVector("nameOfFVector")
  , _nameOfFVectorCompare("nameOfFVectorCompare")
  ,_kMatrDynComp(new KMatrixDynamics(_kMatCompareName, fsParticles, mother2, pathToConfigParserComp, channelId, projectionParticleNames))
{
  addOneGrandMa(_nameOfFVector);
  _FVector = _fVecMap.at(0);
  _kMatrDynComp->addOneGrandMa(_nameOfFVectorCompare);
}

FVectorCompareDynamics::~FVectorCompareDynamics()
{
}

complex<double> FVectorCompareDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  complex<double> result(0.,0.);

  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=_tMatr->kMatrix()->phaseSpaceVec();
  double currentMass=theData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::M_PIPISCAT_NAME));
  _FVector->evalMatrix(currentMass);
  _FVectorCompare->evalMatrix(currentMass);
  evalPhaseCompare(theData, currentMass);
  return result;
}

void FVectorCompareDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  KMatrixDynamics::fillDefaultParams(fitPar);
  _kMatrDynComp->fillDefaultParams(fitPar);
}

void FVectorCompareDynamics::fillParamNameList(){
  _paramNameList.clear();
  KMatrixDynamics::fillParamNameList();
  _kMatrDynComp->fillParamNameList();
}

bool FVectorCompareDynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){
  return true;
}

void FVectorCompareDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  KMatrixDynamics::updateFitParams(fitPar);
  _kMatrDynComp->updateFitParams(fitPar);
}

void FVectorCompareDynamics::fillMasses(EvtData* theData){
  KMatrixDynamics::fillMasses(theData);
  _kMatrDynComp->fillMasses(theData);
}

void FVectorCompareDynamics::evalPhaseCompare(EvtData* theData, double currentMass){
  double currentPhase=std::arg((*_FVector)(_decProjectionIndex, 0));
  double currentPhaseCompare=std::arg((*_FVectorCompare)(_decProjectionIndex, 0));
  double currentPhaseDiff=(currentPhase-currentPhaseCompare)*PawianConstants::radToDeg;
  while(currentPhaseDiff > 180. ) currentPhaseDiff-=180.;
  while(currentPhaseDiff < -180.) currentPhaseDiff+=180.;
  theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME))=currentPhaseDiff;
}


