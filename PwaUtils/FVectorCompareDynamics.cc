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

FVectorCompareDynamics::FVectorCompareDynamics(std::string& name1, std::string& name2, std::vector<Particle*>& fsParticles, Particle* mother1, Particle* mother2, std::string& pathToConfigParser1, std::string& pathToConfigParser2,std::string& pathToFVecCompareConfigFile, std::string dataType, ChannelID channelId, std::string projectionParticleNames) :
AbsDynamics(name1, fsParticles, mother1)
  ,_projectionCompareIndex(0)
  ,_kMatrDyn1(new KMatrixDynamics(name1, fsParticles, mother1, pathToConfigParser1, channelId, projectionParticleNames))
  ,_kMatrDyn2(new KMatrixDynamics(name2, fsParticles, mother2, pathToConfigParser2, channelId, projectionParticleNames))
{
}

FVectorCompareDynamics::~FVectorCompareDynamics()
{
}

complex<double> FVectorCompareDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  complex<double> result(0.,0.);
  return result;
}

void FVectorCompareDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  _kMatrDyn1->fillDefaultParams(fitPar);
  _kMatrDyn2->fillDefaultParams(fitPar);
}

void FVectorCompareDynamics::fillParamNameList(){
  _paramNameList.clear();
  _kMatrDyn1->fillParamNameList();
  _kMatrDyn2->fillParamNameList();
}

bool FVectorCompareDynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){
  return true;
}

void FVectorCompareDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  _kMatrDyn1->updateFitParams(fitPar);
  _kMatrDyn2->updateFitParams(fitPar);
}

void FVectorCompareDynamics::fillMasses(EvtData* theData){
  _kMatrDyn1->fillMasses(theData);
  _kMatrDyn2->fillMasses(theData);
}

void FVectorCompareDynamics::evalPhaseCompare(EvtData* theData, double currentMass){
}


