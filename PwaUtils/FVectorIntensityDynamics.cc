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

// FVectorIntensityDynamics class definition file. -*- C++ -*-
// Copyright 2019 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/FVectorIntensityDynamics.hh"
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
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/PPoleBarrier.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceFactory.hh"
#include "PwaDynamics/KMatrixFunctions.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "Utils/IdStringMapRegistry.hh"
#include "Utils/PawianConstants.hh"
#include "ConfigParser/ParserBase.hh"
#include "ConfigParser/pipiScatteringParser.hh"
#include "qft++Extension/PawianUtils.hh"

FVectorIntensityDynamics::FVectorIntensityDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser,  std::string baseNameFVector, ChannelID channelID, std::string projectionParticleNames) :
  KMatrixDynamics(name, fsParticles, mother, pathToConfigParser, channelID, projectionParticleNames)
  ,_nameOfAmplitudeMag(baseNameFVector+_pVecSuffix+_kMatName+"AmpMag")
  , _nameOfAmplitudePhi(baseNameFVector+_pVecSuffix+_kMatName+"AmpPhi")
  ,_currentAmplitudeVal(complex<double>(0.,0.))

{
  _nameOfFVector= baseNameFVector;
  addOneGrandMa(_nameOfFVector);
  _FVector = fVector(_nameOfFVector);

  ParserBase* parserBase=GlobalEnv::instance()->Channel(channelID)->parser();
  pipiScatteringParser* thePiPiScatteringParser = dynamic_cast<pipiScatteringParser*>(parserBase);

  std::string prodMomParamsStr=thePiPiScatteringParser->prodMomParamsStr();
  std::stringstream prodMomParamsStringStr;
  prodMomParamsStringStr << prodMomParamsStr;
  
  std::string L2prodStr;
  prodMomParamsStringStr >> L2prodStr;
  _L2prod=atof(L2prodStr.c_str());

  std::string s1ProdStr;
  prodMomParamsStringStr >> s1ProdStr;
  _s1Prod=atof(s1ProdStr.c_str());

  std::string s2ProdStr;
  prodMomParamsStringStr >> s2ProdStr;
  _s2Prod=atof(s2ProdStr.c_str());

  if(_s2Prod > _s1Prod){
    Alert << "wrong s values for the production particles: _s2Prod = "
	  << _s2Prod << " > " << "_s1Prod = " << _s1Prod
	  << endmsg;
    Alert << "requirement: _s2Prod <= _s1Prod !!!!" << endmsg;
    exit(1);
  }
}



FVectorIntensityDynamics::~FVectorIntensityDynamics()
{
}

complex<double> FVectorIntensityDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=_tMatr->kMatrix()->phaseSpaceVec();
  double currentMass=theData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::M_PIPISCAT_NAME));

  complex<double> currentFAmp=_FVector->evalProjMatrix( currentMass, _decProjectionIndex, OrbMom);

  double momQ2L=1.;
  if(_L2prod > 0){
    double currentMassS=currentMass*currentMass;
    if(currentMassS > _s1Prod)
      momQ2L=pow(PawianQFT::breakupMomQDefaultFromS(currentMassS, _s1Prod, _s2Prod), _L2prod);
    else momQ2L=pow(PawianQFT::breakupMomQDefaultFromS(_s1Prod, currentMassS, _s2Prod), _L2prod);
  }
  
  double currentResult = norm(_currentAmplitudeVal*currentFAmp*sqrt(thePhpVecs.at(_prodProjectionIndex)->factor(currentMass).real()*currentMass/2. * momQ2L) );
  
  theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME)) = currentResult;
  return currentFAmp;
}

void FVectorIntensityDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  KMatrixDynamics::fillDefaultParams(fitPar);
  fitPar->Add(_nameOfAmplitudeMag, 1. , 0.2);
  fitPar->Add(_nameOfAmplitudePhi, 0. , 0.2);
}

void FVectorIntensityDynamics::fillParamNameList(){
  _paramNameList.clear();
  KMatrixDynamics::fillParamNameList();
  _paramNameList.push_back(_nameOfAmplitudeMag);
  _paramNameList.push_back(_nameOfAmplitudePhi);
}

bool FVectorIntensityDynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){
  return true;
}

void FVectorIntensityDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  KMatrixDynamics::updateFitParams(fitPar);
   _currentAmplitudeVal = std::polar( std::abs(fitPar->Value(_nameOfAmplitudeMag)), 
                               fitPar->Value(_nameOfAmplitudePhi) );
  // _currentAmplitudeVal = std::polar( fitPar->Value(_nameOfAmplitudeMag), 
  //                             fitPar->Value(_nameOfAmplitudePhi) );
}



