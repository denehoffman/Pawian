//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)  
//                 Xiaoshuai Qin (xqin@ep1.rub.de)  
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

//#include <boost/multi_array.hpp>
#include "KMatrixExtract/TMatrixResidueExtr.hh"
#include "KMatrixExtract/TMatrixExtrFcn.hh"
#include "KMatrixExtract/TMatrixExtrFit.hh"
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
#include "ConfigParser/KMatrixParser.hh"
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
#include "FitParams/ParamFactory.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

using namespace ROOT::Minuit2;

TMatrixResidueExtr::TMatrixResidueExtr(std::string pathToConfigParser, std::string pathToFitParams, std::string sheet, std::string pathToSerialzationFile, std::complex<double> energyBorderMin, std::complex<double> energyBorderMax, std::complex<double> energyStartParams): 
  TMatrixExtrBase(pathToConfigParser, pathToFitParams, sheet)
  ,_pathToSerialzationFile(pathToSerialzationFile)
  ,_energyMin(energyBorderMin)
  ,_energyMax(energyBorderMax)
  ,_energyStart(energyStartParams)
  ,_tMatFit(new TMatrixExtrFit(pathToConfigParser, pathToFitParams, sheet, energyBorderMin, energyBorderMax, energyStartParams) )
{
  _phpVecs=_tMatFit->getPhps();
}

TMatrixResidueExtr::~TMatrixResidueExtr()
{
}

bool TMatrixResidueExtr::GetCovMatrix(){
  std::ifstream serializationStream(_pathToSerialzationFile.c_str());

  if(!serializationStream.is_open()){
	Alert << "Could not open serialization file." << endmsg;
	return false;
	//	exit(0);
  }
  _thePwaCovMatrix = std::shared_ptr<PwaCovMatrix>(new PwaCovMatrix);
  boost::archive::text_iarchive boostInputArchive(serializationStream);

  boostInputArchive >> *_thePwaCovMatrix;
  return true;
}

//std::complex<double> TMatrixResidueExtr::CalcMassWidth(){
std::complex<double> TMatrixResidueExtr::CalcMassWidth(std::shared_ptr<AbsPawianParameters> theFitParams){
  _tMatFit->updateTMatDy(theFitParams);
  TMatrixExtrFcn fitFcn(_tMatFit);

  // Set user parameters for MinuitFitFcn
  MnUserParameters upar;
  upar.Add("eReal", _energyStart.real(), 0.001, _energyMin.real(), _energyMax.real());
  upar.Add("eImag", _energyStart.imag(), 0.001, _energyMin.imag(), _energyMax.imag());

  MnMigrad migrad(fitFcn, upar);
  InfoMsg <<"Start Migrad "<< endmsg;
  FunctionMinimum min = migrad();

  // MnMigrad migrad2(fitFcn, min.UserState(), MnStrategy(2));
  // min = migrad2();

  // MnMigrad migrad3(fitFcn, min.UserState(), MnStrategy(3));
  // min = migrad3();

  if(!min.IsValid()) {
	// Try with higher strategy
	InfoMsg <<"FM is invalid, try with strategy = 2."<< endmsg;
	MnMigrad migrad2(fitFcn, min.UserState(), MnStrategy(2));
	min = migrad2();
  }

  // Save final fit parameters and their errors in variables
  double final_eReal = min.UserState().Value("eReal");
  double final_eImag = min.UserState().Value("eImag");

  InfoMsg << "\n\n**************** Minuit FunctionMinimum information ******************" << endmsg;
  if(min.IsValid())             InfoMsg << "\n Function minimum is valid.\n";
  else                          InfoMsg << "\n WARNING: Function minimum is invalid!" << endmsg;
  if(min.HasValidCovariance())  InfoMsg << "\n Covariance matrix is valid." << endmsg;
  else                          InfoMsg << "\n WARNING: Covariance matrix is invalid!" << endmsg;
  InfoMsg <<" # of function calls: " << min.NFcn() << endmsg;
  InfoMsg <<" minimum edm: " << std::setprecision(10) << min.Edm()<<endmsg;
  if(!min.HasValidParameters()) InfoMsg << " hasValidParameters() returned FALSE" << endmsg;
  if(!min.HasAccurateCovar())   InfoMsg << " hasAccurateCovar() returned FALSE" << endmsg;
  if(!min.HasPosDefCovar()){    InfoMsg << " hasPosDefCovar() returned FALSE" << endmsg;
	if(min.HasMadePosDefCovar()) InfoMsg << " hasMadePosDefCovar() returned TRUE" << endmsg;
  }
  if(!min.HasCovariance())      InfoMsg << " hasCovariance() returned FALSE" << endmsg;
  if(min.HasReachedCallLimit()) InfoMsg << " hasReachedCallLimit() returned TRUE" << endmsg;
  if(min.IsAboveMaxEdm())       InfoMsg << " isAboveMaxEdm() returned TRUE" << endmsg;
  if(min.HesseFailed())         InfoMsg << " hesseFailed() returned TRUE\n" << endmsg;
  return std::complex<double>(final_eReal, final_eImag);
}

void TMatrixResidueExtr::Calculation(){

  std::vector<ResidueProperties> resPropAverage;
  std::vector<ResidueProperties> resPropReal;
  std::vector<ResidueProperties> resPropImag;
  std::complex<double> polePos;
  CalcResidueAll(_params, polePos, resPropReal, resPropImag, resPropAverage);

  if (GetCovMatrix() ){
	//  unsigned int nPar = _params->Params().size();
	std::shared_ptr<AbsPawianParameters> newFitParams = std::shared_ptr<AbsPawianParameters>(_params->Clone());
	unsigned int nKMatrixPar = _kMatrixParamNames.size();
	InfoMsg << "KMatrix Params Size " << nKMatrixPar << endmsg;
	double stepSize = 0.0001;
	std::vector<ResidueProperties> resPropStep;
	std::vector<ResidueProperties> resPropTemp;
	std::complex<double> polePosStep;
	std::map< std::string, std::vector<ResidueProperties> > _derivatives;
	for(unsigned int i=0; i<nKMatrixPar; i++){
	  std::string parName = _kMatrixParamNames[i];
	  unsigned int index = _params->Index(parName);
	  double parOrig = _params->Value(index);
	  newFitParams->SetValue(index, parOrig + stepSize);
	  CalcResidueAll(newFitParams, polePosStep, resPropTemp, resPropTemp, resPropStep);
	  std::vector<ResidueProperties> resPropDerivative;
	  resPropDerivative.resize(_phpVecs.size());
	  for (unsigned int i=0; i<resPropStep.size(); ++i){
		resPropDerivative.at(i).absR = (resPropStep.at(i).absR - resPropAverage.at(i).absR)/stepSize;
		resPropDerivative.at(i).theta = (fmod((resPropStep.at(i).theta - resPropAverage.at(i).theta)+5*M_PI, 2*M_PI) - M_PI)/stepSize;
		resPropDerivative.at(i).gammai = (resPropStep.at(i).gammai - resPropAverage.at(i).gammai)/stepSize;
	  }
	  _derivatives.insert(make_pair(parName, resPropDerivative));
	  newFitParams->SetValue(index, parOrig);
	}
	std::vector<ResidueProperties> resPropError;
	resPropError.resize(_phpVecs.size());
	for(unsigned int i=0; i<nKMatrixPar; i++){
	  for(unsigned int j=0; j<nKMatrixPar; j++){
		std::string name1 = _kMatrixParamNames[i];
		std::string name2 = _kMatrixParamNames[j];
		for (unsigned int index=0;index!=resPropAverage.size();index++){
		  resPropError.at(index).absR += (_derivatives[name1].at(index).absR ) * _thePwaCovMatrix->GetElement(name1, name2) * (_derivatives[name2].at(index).absR );
		  resPropError.at(index).theta += (_derivatives[name1].at(index).theta ) * _thePwaCovMatrix->GetElement(name1, name2) * (_derivatives[name2].at(index).theta );
		  resPropError.at(index).gammai += (_derivatives[name1].at(index).gammai ) * _thePwaCovMatrix->GetElement(name1, name2) * (_derivatives[name2].at(index).gammai );
		}
	  }
	}
	for (unsigned int index=0;index!=resPropAverage.size();index++){
	  resPropAverage.at(index).errAbsR = sqrt(resPropError.at(index).absR);
	  resPropAverage.at(index).errTheta = sqrt(resPropError.at(index).theta);
	  resPropAverage.at(index).errGammai = sqrt(resPropError.at(index).gammai);
	}
  }
  dumpResult(polePos, resPropReal, resPropImag, resPropAverage);
  return;
}

void TMatrixResidueExtr::CalcResidueAll(std::shared_ptr<AbsPawianParameters> theFitParams, std::complex<double>& polePos, std::vector<ResidueProperties>& resPropReal, std::vector<ResidueProperties>& resPropImag, std::vector<ResidueProperties>& resPropAverage){
  //  std::vector<ResidueProperties> resPropReal;
  resPropReal.resize(_phpVecs.size());
  //  std::vector<ResidueProperties> resPropImag;
  resPropImag.resize(_phpVecs.size());
  //  vector<ResidueProperties> resPropAverage;
  resPropAverage.resize(_phpVecs.size());

  std::complex<double> result(0.,0.);
  const double epsilon=0.000001;
  polePos = CalcMassWidth(theFitParams);
  std::complex<double> polePosEpsilonImagp = polePos + std::complex<double>(0., epsilon);
  std::complex<double> polePosEpsilonImagm = polePos + std::complex<double>(0., -epsilon);
  std::complex<double> polePosEpsilonRealp = polePos + std::complex<double>(epsilon, 0.);
  std::complex<double> polePosEpsilonRealm = polePos + std::complex<double>(-epsilon, 0.);

  std::shared_ptr<TMatrixRel> currentTMatRealp=_tMatFit->getNewTMat();
  //  currentTMatRealp->evalNonRelMatrix(polePosEpsilonRealp);
  currentTMatRealp->evalMatrix(polePosEpsilonRealp);

  std::shared_ptr<TMatrixRel> currentTMatRealm=_tMatFit->getNewTMat();
  //  currentTMatRealm->evalNonRelMatrix(polePosEpsilonRealm);
  currentTMatRealm->evalMatrix(polePosEpsilonRealm);

  std::shared_ptr<TMatrixRel> currentTMatImagp=_tMatFit->getNewTMat();
  //  currentTMatImagp->evalNonRelMatrix(polePosEpsilonImagp);
  currentTMatImagp->evalMatrix(polePosEpsilonImagp);

  std::shared_ptr<TMatrixRel> currentTMatImagm=_tMatFit->getNewTMat();
  //  currentTMatImagm->evalNonRelMatrix(polePosEpsilonImagm);
  currentTMatImagm->evalMatrix(polePosEpsilonImagm);

  InfoMsg << "\n\nm - i/2. Gamma: " << polePos.real()  << " - i/2. " << -2.*polePos.imag() << endmsg;


  double gammaSum=0.;
  double brSum=0.;  
  for(unsigned int i=0 ; i<_phpVecs.size(); ++i){
	ResidueProperties  currentResPropReal;
	ResidueProperties  currentResPropImag;
	ResidueProperties  currentResPropAverage;

	std::complex<double> resultEpsilonRealpInv= 1./(sqrt(_phpVecs.at(i)->factor(polePos))*(*currentTMatRealp)(i,i)*sqrt(_phpVecs.at(i)->factor(polePos)));
	std::complex<double> resultEpsilonRealmInv= 1./(sqrt(_phpVecs.at(i)->factor(polePos))*(*currentTMatRealm)(i,i)*sqrt(_phpVecs.at(i)->factor(polePos)));
	std::complex<double> resultApproxReal=(resultEpsilonRealpInv-resultEpsilonRealmInv)/(2.*epsilon);
	InfoMsg << "resultApproxRel: " << resultApproxReal << endmsg;

	std::complex<double> resultEpsilonImagpInv= 1./(sqrt(_phpVecs.at(i)->factor(polePos))*(*currentTMatImagp)(i,i)*sqrt(_phpVecs.at(i)->factor(polePos)));
	std::complex<double> resultEpsilonImagmInv= 1./(sqrt(_phpVecs.at(i)->factor(polePos))*(*currentTMatImagm)(i,i)*sqrt(_phpVecs.at(i)->factor(polePos)));    
	std::complex<double> resultApproxImag=1./PawianConstants::i*(resultEpsilonImagpInv-resultEpsilonImagmInv)/(2.*epsilon);
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

	resPropReal.at(i)=currentResPropReal;
	resPropImag.at(i)=currentResPropImag;
	resPropAverage.at(i)=currentResPropAverage;
  }
}

void TMatrixResidueExtr::dumpResult(std::complex<double> polePos, std::vector<ResidueProperties> resPropReal, std::vector<ResidueProperties> resPropImag, std::vector<ResidueProperties> resPropAv){
  std::string oFileName("residues.out");
  std::ofstream theStream(oFileName);
  double gammaTotalr=0.;
  double gammaTotali=0.;
  double gammaTotalav=0.;
  double gammaTotalavError=0.;
  theStream << "m - i/2. Gamma = " << polePos.real()  << " - i/2. " << -2.*polePos.imag() << std::endl << std::endl;
  theStream << setw(7) << "channel" 
	<< setw(12) << "Theta(r)"  << setw(12) << "Gamma_i(r)" << setw(12) << "BR[%](r)"
	<< setw(12) << "Theta(i)"  << setw(12) << "Gamma_i(i)" << setw(12) << "BR[%](i)" 
	<< setw(12) << "Theta(av)" << setw(29) << "Gamma_i(av)" << setw(29) << "BR[%](av)"  << std::endl;

  for (unsigned int i=0; i<resPropReal.size(); ++i){
	theStream << setw(7) << i  
	  << setw(12) << resPropReal.at(i).theta*180./M_PI << setw(12) << resPropReal.at(i).gammai << setw(12) << -resPropReal.at(i).gammai/(2.*polePos.imag())
	  << setw(12) << resPropImag.at(i).theta*180./M_PI << setw(12) << resPropImag.at(i).gammai << setw(12) << -resPropImag.at(i).gammai/(2.*polePos.imag())
	  << setw(12) << resPropAv.at(i).theta*180./M_PI << setw(12) << resPropAv.at(i).gammai << " +/- " << setw(12) << resPropAv.at(i).errGammai 
	  << setw(12) << -resPropAv.at(i).gammai/(2.*polePos.imag()) <<  " +/- " << setw(12) << -resPropAv.at(i).errGammai/(2.*polePos.imag())
	  << std::endl;
	gammaTotalr+=resPropReal.at(i).gammai;
	gammaTotali+=resPropImag.at(i).gammai;
	gammaTotalav+=resPropAv.at(i).gammai; 
	gammaTotalavError+=resPropAv.at(i).errGammai*resPropAv.at(i).errGammai; 
  }  
  theStream << setw(7) << "all" 
	<< setw(36) << -gammaTotalr/(2.*polePos.imag())
	<< setw(36) << -gammaTotali/(2.*polePos.imag())
	<< setw(53) << -gammaTotalav/(2.*polePos.imag()) << " +/- " 
	<< setw(12) << -sqrt(gammaTotalavError)/(2.*polePos.imag()) 
	<< std::endl;
}
