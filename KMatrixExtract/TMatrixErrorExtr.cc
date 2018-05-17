//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)  
//                 Xiaoshuai Qin (xqin@ep1.rub.de)             //
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
#include <vector>
#include <complex>
#include <map>
#include <getopt.h>
#include <iomanip>      // std::setprecision
#include <memory>
#include <algorithm> 

//#include <boost/multi_array.hpp>
#include "KMatrixExtract/TMatrixErrorExtr.hh"
#include "KMatrixExtract/TMatrixExtrFcn.hh"
#include "KMatrixExtract/TMatrixExtrFit.hh"
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

TMatrixErrorExtr::TMatrixErrorExtr(std::string pathToConfigParser, std::string pathToFitParams, std::string sheet, std::string pathToSerialzationFile, std::complex<double> energyBorderMin, std::complex<double> energyBorderMax, std::complex<double> energyStartParams): 
  TMatrixExtrBase(pathToConfigParser, pathToFitParams, sheet)
  ,_pathToSerialzationFile(pathToSerialzationFile)
  ,_energyMin(energyBorderMin)
  ,_energyMax(energyBorderMax)
  ,_energyStart(energyStartParams)
  ,_tMatFit(new TMatrixExtrFit(pathToConfigParser, pathToFitParams, sheet, energyBorderMin, energyBorderMax, energyStartParams) )
{
  //  _tMatFit = new TMatrixExtrFit(pathToConfigParser, pathToFitParams, energyMin, energyMax, energyStart);
}

TMatrixErrorExtr::~TMatrixErrorExtr()
{
}

void TMatrixErrorExtr::GetCovMatrix(){
  std::ifstream serializationStream(_pathToSerialzationFile.c_str());

  if(!serializationStream.is_open()){
	Alert << "Could not open serialization file." << endmsg;
	exit(0);
  }
  _thePwaCovMatrix = std::shared_ptr<PwaCovMatrix>(new PwaCovMatrix);
  boost::archive::text_iarchive boostInputArchive(serializationStream);

  boostInputArchive >> *_thePwaCovMatrix;
}

void TMatrixErrorExtr::CalcOriginal(){
  _result = CalcMassWidth(_params);
  return;
}


void TMatrixErrorExtr::Calculation(){
  GetCovMatrix();
  CalcOriginal();
  std::complex<double> resultErr=0;
  double stepSize = 0.0001;

  //  unsigned int nPar = _params->Params().size();
  std::shared_ptr<AbsPawianParameters> newFitParams = std::shared_ptr<AbsPawianParameters>(_params->Clone());
  unsigned int nKMatrixPar = _kMatrixParamNames.size();
  InfoMsg << "KMatrix Params Size " << nKMatrixPar << endmsg;

  for(unsigned int i=0; i<nKMatrixPar; i++){

	std::string parName = _kMatrixParamNames[i];
	unsigned int index = _params->Index(parName);
	double parOrig = _params->Value(index);

	newFitParams->SetValue(index, parOrig + stepSize);

	std::complex<double> resultDelta = CalcMassWidth(newFitParams);
	std::complex<double> newDerivative = (resultDelta - _result) / stepSize;
	_derivatives[parName] = newDerivative;
	_realDerivatives.push_back(std::make_pair( parName, newDerivative.real() ) );
	_imagDerivatives.push_back(std::make_pair( parName, newDerivative.imag() ) );

	newFitParams->SetValue(index, parOrig);
  }

  for(unsigned int i=0; i<nKMatrixPar; i++){
	for(unsigned int j=0; j<nKMatrixPar; j++){
	  std::string name1 = _kMatrixParamNames[i];
	  std::string name2 = _kMatrixParamNames[j];
	  double err_real = std::real(_derivatives[name1] ) * _thePwaCovMatrix->GetElement(name1, name2) * std::real(_derivatives[name2] );
	  double err_imag = std::imag(_derivatives[name1] ) * _thePwaCovMatrix->GetElement(name1, name2) * std::imag(_derivatives[name2] );
	  std::string name = name1 + " " + name2;
	  _realError.push_back(std::make_pair( name, err_real ) );
	  _imagError.push_back(std::make_pair( name, err_imag ) );
	  resultErr += std::complex<double> (err_real, err_imag );
	}
  }

  resultErr = std::complex<double>(sqrt(std::real(resultErr) ), sqrt(std::imag(resultErr) ) );
  _error = resultErr;
  return;
}


void TMatrixErrorExtr::printErrors(){
  InfoMsg << "TMatrixErrorExtr::printErrors " << endmsg;
  InfoMsg << "Derivatives: " << endmsg;
  std::sort(_realDerivatives.begin(), _realDerivatives.end(), TMatrixErrorExtr::cmp );
  std::sort(_imagDerivatives.begin(), _imagDerivatives.end(), TMatrixErrorExtr::cmp );
  std::sort(_realError.begin(), _realError.end(), TMatrixErrorExtr::cmp );
  std::sort(_imagError.begin(), _imagError.end(), TMatrixErrorExtr::cmp );
  std::vector<std::pair<string, double> >::iterator iter;
  unsigned int shift = _realDerivatives.size() > 5 ? 5 : 1; 
  for (iter=_realDerivatives.begin();iter!=_realDerivatives.begin()+shift;iter++) {
	InfoMsg << "realDerivative " << iter->first << "\t" << iter->second << endmsg;
  }
  for (iter=_imagDerivatives.begin();iter!=_imagDerivatives.begin()+shift;iter++) {
	InfoMsg << "imagDerivative " << iter->first << "\t" << iter->second << endmsg;
  }
  for (iter=_realError.begin();iter!=_realError.begin()+shift;iter++) {
	InfoMsg << "realError " << iter->first << "\t" << iter->second << endmsg;
  }
  for (iter=_imagError.begin();iter!=_imagError.begin()+shift;iter++) {
	InfoMsg << "imagError " << iter->first << "\t" << iter->second << endmsg;
  }
}

std::complex<double> TMatrixErrorExtr::CalcMassWidth(std::shared_ptr<AbsPawianParameters> theFitParams){
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
