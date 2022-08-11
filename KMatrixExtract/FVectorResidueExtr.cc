//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)  
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
//FVectorResidueExtr class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <complex>
#include <map>
#include <getopt.h>
#include <iomanip>      // std::setprecision
#include <memory>
#include <algorithm> 
#include <cmath>

#include "KMatrixExtract/FVectorResidueExtr.hh"
#include "KMatrixExtract/TMatrixExtrFcn.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh" 
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"
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

#include "ConfigParser/KMatrixParser.hh"
#include "ConfigParser/pipiScatteringParser.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Particle/PdtParser.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "FitParams/AbsPawianParamStreamer.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "FitParams/PwaCovMatrix.hh"
#include "Utils/PawianConstants.hh"

#include "ConfigParser/ParserBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/TMatrixDynamics.hh"
#include "PwaUtils/FVectorIntensityDynamics.hh"

#include "FitParams/ParamFactory.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

using namespace ROOT::Minuit2;

FVectorResidueExtr::FVectorResidueExtr(pipiScatteringParser* theParser) :
  TMatrixResidueExtr(theParser)
  ,_pVecName("")
  ,_extractionMethod("Laurent")
{
  init();
}

FVectorResidueExtr::~FVectorResidueExtr()
{
}

void FVectorResidueExtr::init() {

  std::string baseNameFVector=_pipiScatteringParser->baseNameFVector();
  InfoMsg << "baseNameFVector: " << baseNameFVector << endmsg;
  std::vector<std::string> prodSuffixes= _pipiScatteringParser->addKmatrixProdSuffix();
  
  ChannelID channelID(0);
  
   std::string dummyName="dummy";

  _fVectorIntensityDynamics = 
    std::shared_ptr<FVectorIntensityDynamics>(new FVectorIntensityDynamics(dummyName, _fsParticles, _motherParticle, _pathToKMatrixParser, baseNameFVector, channelID, _projectionParticleNames));
  
  _decProjectionIndex = _fVectorIntensityDynamics->decProjectionIndex();
  InfoMsg << "_decProjectionIndex: " << _decProjectionIndex << endmsg;

  _phpVecCurrent = _phpVecs.at(_decProjectionIndex);
  //_fVector=_fVectorIntensityDynamics->getFVector();
  fillParams();
  _fVector=_fVectorIntensityDynamics->getFVector();
  _fVector->SetBumImPartSigns(_signs); 
}


void FVectorResidueExtr::fillParams(){
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
  //if(_pathToFitParams != "") _fVectorIntensityDynamics->updateFitParams(_params);
  _fVecParamNames=_fVectorIntensityDynamics->paramNames();
  for (unsigned int i=0; i<_fVecParamNames.size(); ++i){
    InfoMsg << "F-Vector ParamNames: " << _fVecParamNames.at(i) << endmsg;
  }
  _fVectorIntensityDynamics->fillDefaultParams(params);
  //if(_pathToFitParams != "") _fVectorIntensityDynamics->updateFitParams(_params);
}

void FVectorResidueExtr::updateFMatDy(std::shared_ptr<AbsPawianParameters> params) {
  //updateTMatDy(params);
  _fVectorIntensityDynamics->updateFitParams(params);
  return;
}

void FVectorResidueExtr::Calculation(){
    
  std::vector<ResidueProperties> resPropAverage;
  std::vector<ResidueProperties> resPropReal;
  std::vector<ResidueProperties> resPropImag;
  std::complex<double> polePos;
  CalcResidueAll(_params, polePos, resPropReal, resPropImag, resPropAverage);
  std::vector<std::string> paramNames=_fVecParamNames;
  InfoMsg << "paramNames.size(): " << paramNames.size() << endmsg;
  CalculationError(_fVecParamNames, resPropAverage);

  dumpResult(polePos, resPropReal, resPropImag, resPropAverage);
  printResults(polePos, resPropAverage);
}

void FVectorResidueExtr::CalcResidueAll(std::shared_ptr<AbsPawianParameters> theFitParams, 
					std::complex<double>& polePos, std::vector<ResidueProperties>& resPropReal, 
					std::vector<ResidueProperties>& resPropImag, 
					std::vector<ResidueProperties>& resPropAverage){

  std::vector<ResidueProperties> resPropAverageTMat;
  std::vector<ResidueProperties> resPropRealTMat;
  std::vector<ResidueProperties> resPropImagTMat;
  std::complex<double> polePosTMat;
  TMatrixResidueExtr::CalcResidueAll(theFitParams, polePosTMat, resPropRealTMat, resPropImagTMat, resPropAverageTMat);

  
  //  std::vector<ResidueProperties> resPropReal;
  resPropReal.resize(_phpVecs.size());
  //  std::vector<ResidueProperties> resPropImag;
  resPropImag.resize(_phpVecs.size());
  //  vector<ResidueProperties> resPropAverage;
  resPropAverage.resize(_phpVecs.size());

  std::complex<double> result(0.,0.);
  const double epsilon=0.00001;
  updateFMatDy(theFitParams);
  polePos = CalcMassWidth(theFitParams);

  std::complex<double> polePosEpsilonImagp = polePos + std::complex<double>(0., epsilon);
  std::complex<double> polePosEpsilonImagm = polePos + std::complex<double>(0., -epsilon);
  std::complex<double> polePosEpsilonRealp = polePos + std::complex<double>(epsilon, 0.);
  std::complex<double> polePosEpsilonRealm = polePos + std::complex<double>(-epsilon, 0.);


  std::vector< complex<double> > currentFVecRealp;
  currentFVecRealp.resize(_phpVecs.size()); 
  _fVector->evalMatrix(polePosEpsilonRealp, _orbitalL);
  for(unsigned int i=0 ; i<_phpVecs.size(); ++i) {
    currentFVecRealp.at(i)= (*_fVector)(i,0);
  }

  std::vector< complex<double> > currentFVecRealm;
  currentFVecRealm.resize(_phpVecs.size()); 
  _fVector->evalMatrix(polePosEpsilonRealm, _orbitalL);
  for(unsigned int i=0 ; i<_phpVecs.size(); ++i) {
    currentFVecRealm.at(i)= (*_fVector)(i,0);
  }


  std::vector< complex<double> > currentFVecImagp;
  currentFVecImagp.resize(_phpVecs.size()); 
  _fVector->evalMatrix(polePosEpsilonImagp, _orbitalL);
  for(unsigned int i=0 ; i<_phpVecs.size(); ++i) {
    currentFVecImagp.at(i)= (*_fVector)(i,0);
  }

  std::vector< complex<double> > currentFVecImagm;
  currentFVecImagm.resize(_phpVecs.size()); 
  _fVector->evalMatrix(polePosEpsilonImagm, _orbitalL);
  for(unsigned int i=0 ; i<_phpVecs.size(); ++i) {
    currentFVecImagm.at(i)= (*_fVector)(i,0);
  }


  InfoMsg << "\n\nm - i/2. Gamma: " << polePos.real()  << " - i/2. " << -2.*polePos.imag() << endmsg;

  for(unsigned int i=0 ; i<_phpVecs.size(); ++i) {
    ResidueProperties  currentResPropReal;
    ResidueProperties  currentResPropImag;
    ResidueProperties  currentResPropAverage;
    
    std::complex<double> resultEpsilonRealpInv= 1./(currentFVecRealp.at(i)*sqrt(_phpVecs.at(i)->factor(polePos)));
    std::complex<double> resultEpsilonRealmInv= 1./(currentFVecRealm.at(i)*sqrt(_phpVecs.at(i)->factor(polePos)));
    std::complex<double> resultApproxReal=(resultEpsilonRealpInv-resultEpsilonRealmInv)/(2.*epsilon);
    InfoMsg << "currentFVecRealp.at(" << i << "): " << currentFVecRealp.at(i) << endmsg;
    InfoMsg << "currentFVecRealm.at(" << i << "): " << currentFVecRealm.at(i) << endmsg;
    InfoMsg << "_phpVecs.at(i)->factor(polePos): " << _phpVecs.at(i)->factor(polePos) << endmsg; 
    InfoMsg << "resultApproxRel: " << resultApproxReal << endmsg;
    
    std::complex<double> resultEpsilonImagpInv= 1./(currentFVecImagp.at(i)*sqrt(_phpVecs.at(i)->factor(polePos)));
    std::complex<double> resultEpsilonImagmInv= 1./(currentFVecImagm.at(i)*sqrt(_phpVecs.at(i)->factor(polePos)));
    std::complex<double> resultApproxImag=1./PawianConstants::i *
      (resultEpsilonImagpInv-resultEpsilonImagmInv)/(2.*epsilon);

    InfoMsg << "currentFVecImagp.at(" << i << "): " << currentFVecImagp.at(i) << endmsg;
    InfoMsg << "currentFVecImagm.at(" << i << "): " << currentFVecImagm.at(i) << endmsg;
    InfoMsg << "resultApproxImag: " << resultApproxImag << endmsg;
    
    std::complex<double> resultApprox = (resultApproxReal+resultApproxImag)/2.;
    InfoMsg << "resultApprox: " << resultApprox << endmsg;
    
    currentResPropReal.absR=abs(1./resultApproxReal);
    currentResPropImag.absR=abs(1./resultApproxImag);
    currentResPropAverage.absR=abs(1./resultApprox);
    
    currentResPropReal.theta=atan2(imag(1./resultApproxReal),real(1./resultApproxReal));
    currentResPropImag.theta=atan2(imag(1./resultApproxImag),real(1./resultApproxImag));
    currentResPropAverage.theta=atan2(imag(1./resultApprox),real(1./resultApprox));
    
    currentResPropReal.gammai=2.*abs(1./resultApproxReal);
    currentResPropImag.gammai=2.*abs(1./resultApproxImag);
    currentResPropAverage.gammai=2.*abs(1./resultApprox);    

    currentResPropReal.gammaiBRj=2.*abs(1./resultApproxReal)*2.*abs(1./resultApproxReal)/(-2.*polePos.imag());
    currentResPropImag.gammaiBRj=2.*abs(1./resultApproxImag)*2.*abs(1./resultApproxImag)/(-2.*polePos.imag());
    currentResPropAverage.gammaiBRj=2.*abs(1./resultApprox)*2.*abs(1./resultApprox)/(-2.*polePos.imag());

    currentResPropReal.gammaProd=2.*abs(1./resultApproxReal)*2.*abs(1./resultApproxReal)/resPropAverageTMat.at(i).gammai;
    currentResPropImag.gammaProd=2.*abs(1./resultApproxImag)*2.*abs(1./resultApproxImag)/resPropAverageTMat.at(i).gammai;
    currentResPropAverage.gammaProd=2.*abs(1./resultApprox)*2.*abs(1./resultApprox)/resPropAverageTMat.at(i).gammai;

    
    resPropReal.at(i)=currentResPropReal;
    resPropImag.at(i)=currentResPropImag;
    resPropAverage.at(i)=currentResPropAverage;
  }
}

void FVectorResidueExtr::dumpResult(std::complex<double> polePos, std::vector<ResidueProperties> resPropReal,
                                    std::vector<ResidueProperties> resPropImag,
                                    std::vector<ResidueProperties> resPropAv){
  std::string oFileName="residuesFVector_"+_extractionMethod+".out";
  std::ofstream theStream(oFileName);
  double gammaTotalr=0.;
  double gammaTotali=0.;
  double gammaTotalav=0.;
  double gammaTotalavError=0.;
  theStream << "m - i/2. Gamma = " << polePos.real()  << " - i/2. " << -2.*polePos.imag() << std::endl << std::endl;
  for (unsigned int i=0; i<resPropReal.size(); ++i){
    theStream << "\nchannel :              " << i << std::endl
              << "\nTheta_r :              "  << resPropReal.at(i).theta*180./M_PI 
	      << "\nTheta_i :              "  << resPropImag.at(i).theta*180./M_PI 
	      << "\nTheta_av:              "  << resPropAv.at(i).theta*180./M_PI 
              << "\nGammaProd*Gamma_i(r) : " << resPropReal.at(i).gammai*resPropReal.at(i).gammai
              << "\nGammaProd*Gamma_i(i) : " << resPropImag.at(i).gammai*resPropImag.at(i).gammai
      	      << "\nGammaProd*Gamma_i(av): " << resPropAv.at(i).gammai*resPropAv.at(i).gammai
	      << " +/- " << 2.*resPropAv.at(i).gammai*resPropAv.at(i).errGammai
	      << "\nGammaProd*BR[%]_i(r)   : " << resPropReal.at(i).gammaiBRj
	      << "\nGammaProd*BR[%]_i(i)   : " << resPropImag.at(i).gammaiBRj
	      << "\nGammaProd*BR[%]_i(av)   : " << resPropAv.at(i).gammaiBRj
	      <<  " +/- " << resPropAv.at(i).errGammaiBRj
	      <<"\nGammaProd(r)             : " << resPropReal.at(i).gammaProd
              <<"\nGammaProd(i)             : "	<< resPropImag.at(i).gammaProd
      	      <<"\nGammaProd(av)            : " << resPropAv.at(i).gammaProd
	      <<  " +/- " << resPropAv.at(i).errGammaProd      
              << std::endl;

        gammaTotalr+=resPropReal.at(i).gammai*resPropReal.at(i).gammai;
    gammaTotali+=resPropImag.at(i).gammai*resPropImag.at(i).gammai;
    gammaTotalav+=resPropAv.at(i).gammai*resPropAv.at(i).gammai;
  }
    theStream << "\nall"
            << "\ngammaTotalr: " << gammaTotalr
            << "\ngammaTotalr: " << gammaTotali
            << "\ngammaTotalav: " << gammaTotalav << " +/- "
            << "\ngammaTotalError: " << gammaTotalavError
            << std::endl;
}

//void FVectorResidueExtr::dumpResult(std::complex<double> polePos, std::vector<ResidueProperties> resPropReal, 
//				    std::vector<ResidueProperties> resPropImag, 
//				    std::vector<ResidueProperties> resPropAv){
//  std::string oFileName="residuesFVector_"+_extractionMethod+".out";
//  std::ofstream theStream(oFileName);
//  double gammaTotalr=0.;
//  double gammaTotali=0.;
//  double gammaTotalav=0.;
//  double gammaTotalavError=0.;
//  theStream << "m - i/2. Gamma = " << polePos.real()  << " - i/2. " << -2.*polePos.imag() << std::endl << std::endl;
//  theStream << setw(7) << "channel" 
//	    << setw(10) << "Theta(r)"
//	    << setw(22) << "GammaProd*Gamma_i(r)"
//	    << setw(22) << "GammaProd*BR[%](r)"
//	    << setw(10) << "Theta(i)"
//	    << setw(22) << "GammaProd*Gamma_i(i)"
//	    << setw(22) << "GammaProd*BR[%](i)" 
//	    << setw(10) << "Theta(av)"
//	    << setw(22) << "GammaProd*Gamma_i(av)"
//	    << setw(22) << "GammaProd*BR[%](av)"
//	    << std::endl;
  
//  for (unsigned int i=0; i<resPropReal.size(); ++i){
//    theStream << setw(7) << i  
//	      << setw(10) << resPropReal.at(i).theta*180./M_PI
//	      << setw(22) << resPropReal.at(i).gammai*resPropReal.at(i).gammai 
//	      << setw(22) << (resPropReal.at(i).gammai*resPropReal.at(i).gammai)/(-2.*polePos.imag())
//	      << setw(10) << resPropImag.at(i).theta*180./M_PI
//	      << setw(22) << resPropImag.at(i).gammai*resPropImag.at(i).gammai 
//	      << setw(22) << (resPropImag.at(i).gammai*resPropImag.at(i).gammai )/(-2.*polePos.imag())
//	      << setw(10) << resPropAv.at(i).theta*180./M_PI
//	      << setw(22) << resPropAv.at(i).gammai*resPropAv.at(i).gammai 
//	      << " +/- "
//	      << setw(8) << 2.*resPropAv.at(i).gammai*resPropAv.at(i).errGammai 
//              << setw(22) << resPropAv.at(i).gammaiBRj <<  " +/- "                              
//	      << setw(8) << resPropAv.at(i).errGammaiBRj
//                    << setw(22) << resPropAv.at(i).gammaProd <<  " +/- "
//              << setw(8) << resPropAv.at(i).errGammaProd
//	      << std::endl;
//    gammaTotalr+=resPropReal.at(i).gammai*resPropReal.at(i).gammai;
//    gammaTotali+=resPropImag.at(i).gammai*resPropImag.at(i).gammai;
//    gammaTotalav+=resPropAv.at(i).gammai*resPropAv.at(i).gammai;
//  }

//  theStream << setw(7) << "all" 
//	    << setw(36) << gammaTotalr
//	    << setw(36) << gammaTotali
//	    << setw(53) << gammaTotalav << " +/- " 
//	    << setw(12) << gammaTotalavError 
//	    << std::endl;
//}
