//************************************************************************//
//									  //
//  Copyright 2025 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// TMatrixSigmaParDynamics class definition file. -*- C++ -*-
// Copyright 2025 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "Utils/PawianIOUtils.hh"
#include "PwaUtils/TMatrixSigmaParDynamics.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataScatteringList.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceFactory.hh"
#include "PwaDynamics/KMatrixFunctions.hh"
#include "PwaDynamics/SigmaParameterization.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "Utils/IdStringMapRegistry.hh"
#include "Utils/PawianConstants.hh"

TMatrixSigmaParDynamics::TMatrixSigmaParDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser, std::string dataType, std::string projectionParticleNames) :
  TMatrixDynamics(name,fsParticles, mother)
  ,_sigmaParameterization(std::shared_ptr<SigmaParameterization>(new SigmaParameterization()))
  ,_currentPoleMass(_sigmaParameterization->barePoleMass())
  ,_currentb1(_sigmaParameterization->b1Param())
  ,_currentb2(_sigmaParameterization->b2Param())
  ,_currenta(_sigmaParameterization->aParam())
  ,_currentg4pi(_sigmaParameterization->g4piParam())
  ,_poleMassName("sigmaPoleBareMass")
  ,_b1Name("sigmab1")
  ,_b2Name("sigmab2")
  ,_aName("sigmaa")
  ,_g4piName("sigmag4pi")
{
  
  if(dataType=="Elasticity") _dataTypeID=1;
  else if(dataType=="Phase") _dataTypeID=2;
  else if(dataType=="ArgandUnits") _dataTypeID=3;
  else if(dataType=="Data") _dataTypeID=4;
  else if(dataType=="PhaseDiff") _dataTypeID=5; //in TMatrixCompareDynamics and FVectorCompareDynamics only
  else if(dataType=="PVecIntensity") _dataTypeID=6; //in FVectorIntensityDynamics only
  else if(dataType=="Treal") _dataTypeID=7;
  else if(dataType=="Timag") _dataTypeID=8;
  else if(dataType=="Tcheck") _dataTypeID=9;
  else if(dataType=="Tabs") _dataTypeID=10;
  else if(dataType=="PhasePhi") _dataTypeID=11;
  else {
    Alert << "production formalism/data type with the name" << dataType 
	  << " is not supported for pi pi scattering fits! \n It is working for: "
	  << "Elasticity, ArgandUnits, Phase or PhaseDiff!!!" << endmsg;
    exit(1); 
  }
}


TMatrixSigmaParDynamics::TMatrixSigmaParDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  TMatrixDynamics(name,fsParticles, mother)
  ,_sigmaParameterization(std::shared_ptr<SigmaParameterization>(new SigmaParameterization()))
  ,_currentPoleMass(_sigmaParameterization->barePoleMass())
  ,_currentb1(_sigmaParameterization->b1Param())
  ,_currentb2(_sigmaParameterization->b2Param())
  ,_currenta(_sigmaParameterization->aParam())
  ,_currentg4pi(_sigmaParameterization->g4piParam())
  ,_poleMassName("sigmaPoleBareMass")
  ,_b1Name("sigmab1")
  ,_b2Name("sigmab2")
  ,_aName("sigmaa")
  ,_g4piName("sigmag4pi")
{
}

TMatrixSigmaParDynamics::~TMatrixSigmaParDynamics()
{
}

complex<double> TMatrixSigmaParDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

  double currentMass =
    theData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::M_PIPISCAT_NAME));

  if(_dataTypeID==1) evalElasticity(theData, currentMass, OrbMom);   
  else if(_dataTypeID==2){
    // if(!_prodIsNotDecChannel) evalPhase(theData, currentMass, OrbMom);
    //    else evalRelativePhase(theData, currentMass, OrbMom);
    evalPhase(theData, currentMass, OrbMom);
  }
  else if(_dataTypeID==3){
    evalArgandUnits(theData, currentMass, OrbMom);
  }
  else if(_dataTypeID==7){
    evalTreal(theData, currentMass, OrbMom);
  }
  else if(_dataTypeID==8){
    evalTimag(theData, currentMass, OrbMom);
  }
  else if(_dataTypeID==10){
    evalTabs(theData, currentMass, OrbMom);
  }
  else if(_dataTypeID==11){
    evalPhasePhi(theData, currentMass, OrbMom);
  }
  else{
    Alert << "_dataTypeID = " <<_dataTypeID << " is not supported!!!" << endmsg;
    exit(1); 
  }

  return complex<double>(1.,0);
}

void TMatrixSigmaParDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  //fill pole bare mass
  fitPar->Add(_poleMassName, _currentPoleMass, 0.01, std::abs(_currentPoleMass-0.2), _currentPoleMass+0.2);
  //fill b1 param
  fitPar->Add(_b1Name, _currentb1, 0.01, std::abs(_currentb1-0.2), _currentb1+0.2);
  //fill b2 param
  fitPar->Add(_b2Name, _currentb2, 0.01, std::abs(_currentb2-0.2), _currentb2+0.2);
  //fill a param
  fitPar->Add(_aName, _currenta, 0.01, std::abs(_currenta-0.2), _currenta+0.2);
  //fill g4pi param
  fitPar->Add(_g4piName, _currentg4pi, 0.00001, 0., _currentg4pi+0.0001);
}

void TMatrixSigmaParDynamics::fillParamNameList(){
  _paramNameList.clear();
  //fill
  _paramNameList.push_back(_poleMassName);
  _paramNameList.push_back(_b1Name);
  _paramNameList.push_back(_b2Name);
  _paramNameList.push_back(_aName);
  _paramNameList.push_back(_g4piName);
}

// bool TMatrixSigmaParDynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){
//   return true;
// }

void TMatrixSigmaParDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){

  _sigmaParameterization->setBarePoleMass(_currentPoleMass);

  _currentb1= fitPar->Value(_b1Name);
  _sigmaParameterization->setb1Param(_currentb1);

  _currentb2= fitPar->Value(_b2Name);
  _sigmaParameterization->setb2Param(_currentb2);

  _currenta= fitPar->Value(_aName);
  _sigmaParameterization->setaParam(_currenta);

  _currentg4pi= fitPar->Value(_g4piName);
  _sigmaParameterization->setg4piParam(_currentg4pi);
}



void TMatrixSigmaParDynamics::evalElasticity(EvtData* theData, double currentMass, Spin OrbMom){
  complex<double> currentTijRel= _sigmaParameterization->calcT(currentMass);
  std::complex<double> rho2pi=_sigmaParameterization->rho2pi();
  
  complex<double> SijRel=complex<double>(1.,0.)+2.*PawianConstants::i *
    sqrt(rho2pi.real()*rho2pi.real()) * currentTijRel;
  
  theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME))=sqrt(norm(SijRel));
}


void TMatrixSigmaParDynamics::evalPhase(EvtData* theData, double currentMass, Spin OrbMom=0){
  complex<double> currentTijRel= _sigmaParameterization->calcT(currentMass);
  std::complex<double> rho2pi=_sigmaParameterization->rho2pi();

  complex<double> currentTiiRel_rho= currentTijRel*rho2pi.real();
  double currentReERel = currentTiiRel_rho.real();
  double currentImERel = currentTiiRel_rho.imag() - 0.5;

  double deltaRel = 0.5*atan2(currentImERel, fabs(currentReERel))*PawianConstants::radToDeg + 45.0;
  if (currentTiiRel_rho.real()  < 0.0) {deltaRel = 180.0 - deltaRel;}

  while(deltaRel>180.) deltaRel-=180.;
  while(deltaRel<0.) deltaRel+=180.;

  double phiData=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::DATA_PIPISCAT_NAME));

  while( (phiData-deltaRel) > 90.) deltaRel+=180.;
  while( (deltaRel-phiData) > 90.) deltaRel-=180.;

  theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME))=deltaRel;

}

void TMatrixSigmaParDynamics::evalRelativePhase(EvtData* theData, double currentMass, Spin OrbMom){
  Alert << "mode evalRelativePhase is not supported!!!" << endmsg;
  exit(1);
}

void TMatrixSigmaParDynamics::evalArgandUnits(EvtData* theData, double currentMass, Spin OrbMom){
  Alert << "mode evalArgandUnits is not supported!!!" << endmsg;
  exit(1);
}

void TMatrixSigmaParDynamics::evalTreal(EvtData* theData, double currentMass, Spin OrbMom){
  Alert << "mode evalTreal is not supported!!!" << endmsg;
  exit(1);
}

void TMatrixSigmaParDynamics::evalTimag(EvtData* theData, double currentMass, Spin OrbMom){
  Alert << "mode evalTimag is not supported!!!" << endmsg;
  exit(1);
}

void TMatrixSigmaParDynamics::evalTabs(EvtData* theData, double currentMass, Spin OrbMom){
  Alert << "mode evalTabs is not supported!!!" << endmsg;
  exit(1);
}

void TMatrixSigmaParDynamics::evalPhasePhi(EvtData* theData, double currentMass, Spin OrbMom=0){
  Alert << "mode evalPhasePhi is not supported!!!" << endmsg;
  exit(1);
}

