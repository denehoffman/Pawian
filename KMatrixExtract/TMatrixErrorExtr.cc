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

#include "KMatrixExtract/TMatrixErrorExtr.hh"
#include "KMatrixExtract/TMatrixExtrFcn.hh"
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
#include "Utils/PawianIOUtils.hh"

#include "ConfigParser/ParserBase.hh"
#include "ConfigParser/pipiScatteringParser.hh"

#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/TMatrixDynamics.hh"
#include "FitParams/ParamFactory.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

using namespace ROOT::Minuit2;

TMatrixErrorExtr::TMatrixErrorExtr(pipiScatteringParser* theParser) :
  TMatrixExtrBase(theParser)
  ,_pathToSerialzationFile(theParser->serializationFile())
  ,_energyMin(complex<double>(theParser->minRealMass(), -std::abs(theParser->minImagMass())))
  ,_energyMax(complex<double>(theParser->maxRealMass(), -std::abs(theParser->maxImagMass())))
  ,_energyStart(complex<double>(theParser->startRealMass(), -std::abs(theParser->startImagMass())))
  ,_calcWithErrors(true)
{
  if(_isFixedKMatrixParametrization){
    _pathToSerialzationFile= PawianIOUtils::getFileName(GlobalEnv::instance()->KMatrixStorePath(), _kMatrixParser->serialFixedParamFile());
  }
}

TMatrixErrorExtr::~TMatrixErrorExtr()
{
}

bool TMatrixErrorExtr::GetCovMatrix(){
  bool result=true;
  std::ifstream serializationStream(_pathToSerialzationFile.c_str());

  if(!serializationStream.is_open()){
	InfoMsg << "Could not open serialization file." << endmsg;
	InfoMsg << "calculation without errors!!!" << endmsg;       
	result=false;
  }
  else{
    _thePwaCovMatrix = std::shared_ptr<PwaCovMatrix>(new PwaCovMatrix);
    boost::archive::text_iarchive boostInputArchive(serializationStream);
    boostInputArchive >> *_thePwaCovMatrix;
  }
  return result;
}

void TMatrixErrorExtr::CalcOriginal(){
  _result = CalcMassWidth(_params);
  return;
}

void TMatrixErrorExtr::Calculation(){
  _calcWithErrors = GetCovMatrix();
  if (_calcWithErrors) CalcWithErrrors();
  else  CalcOriginal();
  
}

void TMatrixErrorExtr::CalcWithErrrors(){
  GetCovMatrix();
  CalcOriginal();
  std::complex<double> resultErr=0;
  double stepSize = 0.0001;

  //  unsigned int nPar = _params->Params().size();
  std::shared_ptr<AbsPawianParameters> newFitParams = std::shared_ptr<AbsPawianParameters>(_params->Clone());
  unsigned int nKMatrixPar = _kMatrixParamNames.size();

  for(unsigned int i=0; i<nKMatrixPar; i++) {

    std::string parName = _kMatrixParamNames[i];
    unsigned int index = _params->Index(parName);
    double parOrig = _params->Value(index);

    stepSize=sqrt(_thePwaCovMatrix->GetElement(parName, parName))/500.;
    if(stepSize<1.e-10) stepSize=1.e-10;
    
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
      double err_real = std::real(_derivatives[name1] ) * 
	_thePwaCovMatrix->GetElement(name1, name2) * std::real(_derivatives[name2] );
      double err_imag = std::imag(_derivatives[name1] ) * 
	_thePwaCovMatrix->GetElement(name1, name2) * std::imag(_derivatives[name2] );
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
  if(!_calcWithErrors) return;
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

std::complex<double> TMatrixErrorExtr::CalcMassWidth(std::shared_ptr<AbsPawianParameters> theFitParams) {
  std::shared_ptr<TMatrixErrorExtr> tmpTMatrixErrorExtr(new TMatrixErrorExtr(_pipiScatteringParser));
  tmpTMatrixErrorExtr->updateTMatDy(theFitParams);
  TMatrixExtrFcn fitFcn(tmpTMatrixErrorExtr);

  // updateTMatDy(theFitParams);
  // TMatrixExtrFcn fitFcn(shared_from_this());

  // Set user parameters for MinuitFitFcn
  MnUserParameters upar;
  InfoMsg << "_energyStart.imag(): " << _energyStart.imag() << endmsg;
  InfoMsg << "_energyMin.imag(): " << _energyMin.imag() << endmsg;
  InfoMsg << "_energyMax.imag(): " << _energyMax.imag() << endmsg;
  upar.Add("eReal", _energyStart.real(), 0.001, _energyMin.real(), _energyMax.real());
  upar.Add("eImag", _energyStart.imag(), 0.001, _energyMin.imag(), _energyMax.imag());

  MnMigrad migrad(fitFcn, upar);
  InfoMsg <<"Start Migrad "<< endmsg;
  FunctionMinimum min = migrad();

  if(!min.IsValid()) {
    // Try with higher strategy
    InfoMsg <<"FM is invalid, try with strategy = 2."<< endmsg;
    MnMigrad migrad2(fitFcn, min.UserState(), MnStrategy(2));
    min = migrad2();
  }

  //start second iteration
  // MnMigrad migrad1a(fitFcn, min.UserState(), MnStrategy(1));
  // min = migrad1a();
  // if(!min.IsValid()) {
  //   // Try with higher strategy
  //   InfoMsg <<"FM is invalid, try with strategy = 2."<< endmsg;
  //   MnMigrad migrad2a(fitFcn, min.UserState(), MnStrategy(2));
  //   min = migrad2a();
  // }
  
  // Save final fit parameters and their errors in variables
  double final_eReal = min.UserState().Value("eReal");
  double final_eImag = min.UserState().Value("eImag");

  InfoMsg << "\n\n**************** Minuit FunctionMinimum information ******************" << endmsg;
  if(min.IsValid()) {
    InfoMsg << "\n Function minimum is valid.\n" << endmsg;
  } else {
    InfoMsg << "\n WARNING: Function minimum is invalid!" << endmsg;
  }
  if(min.HasValidCovariance()) {
    InfoMsg << "\n Covariance matrix is valid." << endmsg;
  } else {
    InfoMsg << "\n WARNING: Covariance matrix is invalid!" << endmsg;
  }
  InfoMsg <<" # of function calls: " << min.NFcn() << endmsg;
  InfoMsg <<" minimum edm: " << std::setprecision(10) << min.Edm() << endmsg;
  if(!min.HasValidParameters()) {
    InfoMsg << " hasValidParameters() returned FALSE" << endmsg;
  }
  if(!min.HasAccurateCovar()) {
    InfoMsg << " hasAccurateCovar() returned FALSE" << endmsg;
  }
  if(!min.HasPosDefCovar()) {
    InfoMsg << " hasPosDefCovar() returned FALSE" << endmsg;
    if(min.HasMadePosDefCovar()) {
      InfoMsg << " hasMadePosDefCovar() returned TRUE" << endmsg;
    }
  }
  if(!min.HasCovariance()) {
    InfoMsg << " hasCovariance() returned FALSE" << endmsg;
  }
  if(min.HasReachedCallLimit()) {
    InfoMsg << " hasReachedCallLimit() returned TRUE" << endmsg;
  }
  if(min.IsAboveMaxEdm()) {
    InfoMsg << " isAboveMaxEdm() returned TRUE" << endmsg;
  }
  if(min.HesseFailed()) {
    InfoMsg << " hesseFailed() returned TRUE\n" << endmsg;
  }

  double final_eRealError = min.UserState().Error("eReal");
  double final_eImagError = min.UserState().Error("eImag");

  InfoMsg << "\n\n****** fit error for pole position  ********" << endmsg;  
  InfoMsg << "Real error: " <<  final_eRealError << endmsg;
  InfoMsg << "eImag error: " <<  final_eImagError << endmsg;
  InfoMsg << "\n***********************\n"  << endmsg;
  return std::complex<double>(final_eReal, final_eImag);
}
