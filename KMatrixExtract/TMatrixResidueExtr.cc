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

#include "KMatrixExtract/TMatrixResidueExtr.hh"
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
#include "FitParams/ParamFactory.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

using namespace ROOT::Minuit2;

TMatrixResidueExtr::TMatrixResidueExtr(pipiScatteringParser* theParser) :
  TMatrixErrorExtr(theParser)
  ,_extractionMethod("Laurent")
{
}

TMatrixResidueExtr::~TMatrixResidueExtr()
{
}

void TMatrixResidueExtr::Calculation(){

  std::vector<ResidueProperties> resPropAverage;
  std::vector<ResidueProperties> resPropReal;
  std::vector<ResidueProperties> resPropImag;
  std::complex<double> polePos;
  CalcResidueAll(_params, polePos, resPropReal, resPropImag, resPropAverage);

  std::vector<std::string> paramNames=_kMatrixParamNames;
  
  CalculationError(paramNames, resPropAverage);

  dumpResult(polePos, resPropReal, resPropImag, resPropAverage);
  printResults(polePos, resPropAverage);
  return;
}

void TMatrixResidueExtr::CalculationError(std::vector<std::string> paramNames, std::vector<ResidueProperties>& resPropAverage){
  if (GetCovMatrix() ) {
    std::shared_ptr<AbsPawianParameters> newFitParams = std::shared_ptr<AbsPawianParameters>(_params->Clone());
    unsigned int nParams = paramNames.size();
    double stepSize = 0.0001;
    std::vector<ResidueProperties> resPropStep;
    std::vector<ResidueProperties> resPropTemp;
    std::complex<double> polePosStep;
    std::map< std::string, std::vector<ResidueProperties> > derivatives;

    for(unsigned int i=0; i<nParams; i++){
      std::string parName = paramNames.at(i);
      unsigned int index = _params->Index(parName);
      double parOrig = _params->Value(index);
      stepSize=sqrt(_thePwaCovMatrix->GetElement(parName, parName))/500.;
      if (stepSize<1.e-10) stepSize=1.e-10;
      newFitParams->SetValue(index, parOrig + stepSize);
      CalcResidueAll(newFitParams, polePosStep, resPropTemp, resPropTemp, resPropStep);
      std::vector<ResidueProperties> resPropDerivative;
      resPropDerivative.resize(_phpVecs.size());
      for (unsigned int i=0; i<resPropStep.size(); ++i){
        resPropDerivative.at(i).absR = (resPropStep.at(i).absR - resPropAverage.at(i).absR)/stepSize;
        resPropDerivative.at(i).theta = (fmod((resPropStep.at(i).theta - 
                                               resPropAverage.at(i).theta)+5*M_PI, 2*M_PI) - M_PI)/stepSize;
        resPropDerivative.at(i).gammai = (resPropStep.at(i).gammai - resPropAverage.at(i).gammai)/stepSize;
	resPropDerivative.at(i).gammaigammaj = (resPropStep.at(i).gammaigammaj - resPropAverage.at(i).gammaigammaj)/stepSize;
	resPropDerivative.at(i).gammaiBRj = (resPropStep.at(i).gammaiBRj - resPropAverage.at(i).gammaiBRj)/stepSize;
      }
      derivatives.insert(make_pair(parName, resPropDerivative));
      newFitParams->SetValue(index, parOrig);
    }

    std::vector<ResidueProperties> resPropError;
    resPropError.resize(_phpVecs.size());
    for(unsigned int i=0; i<nParams; i++){
      for(unsigned int j=0; j<nParams; j++){
        std::string name1 = paramNames.at(i);
        std::string name2 = paramNames.at(j);
        for (unsigned int index=0;index!=resPropAverage.size();index++){
          resPropError.at(index).absR += (derivatives[name1].at(index).absR ) * 
            _thePwaCovMatrix->GetElement(name1, name2) * (derivatives[name2].at(index).absR );
          resPropError.at(index).theta += (derivatives[name1].at(index).theta ) * 
            _thePwaCovMatrix->GetElement(name1, name2) * (derivatives[name2].at(index).theta );
          resPropError.at(index).gammai += (derivatives[name1].at(index).gammai ) * 
            _thePwaCovMatrix->GetElement(name1, name2) * (derivatives[name2].at(index).gammai );
	  resPropError.at(index).gammaigammaj += (derivatives[name1].at(index).gammaigammaj ) * 
            _thePwaCovMatrix->GetElement(name1, name2) * (derivatives[name2].at(index).gammaigammaj );
	  resPropError.at(index).gammaiBRj += (derivatives[name1].at(index).gammai ) *
            _thePwaCovMatrix->GetElement(name1, name2) * (derivatives[name2].at(index).gammaiBRj );
        }
      }
    }
    for (unsigned int index=0;index!=resPropAverage.size();index++){
      resPropAverage.at(index).errAbsR = sqrt(resPropError.at(index).absR);
      resPropAverage.at(index).errTheta = sqrt(resPropError.at(index).theta);
      resPropAverage.at(index).errGammai = sqrt(resPropError.at(index).gammai);
      resPropAverage.at(index).errGammaiGammaj = sqrt(resPropError.at(index).gammaigammaj);
      resPropAverage.at(index).errGammaiBRj = sqrt(resPropError.at(index).gammaiBRj);
    }
  }
}

void TMatrixResidueExtr::CalcResidueAll(std::shared_ptr<AbsPawianParameters> theFitParams, 
					std::complex<double>& polePos, std::vector<ResidueProperties>& resPropReal, 
					std::vector<ResidueProperties>& resPropImag, 
					std::vector<ResidueProperties>& resPropAverage){
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
 
  std::shared_ptr<TMatrixRel> currentTMatRealp = getNewTMat();
  //  currentTMatRealp->evalNonRelMatrix(polePosEpsilonRealp);
  currentTMatRealp->evalMatrix(polePosEpsilonRealp, _orbitalL);

  std::shared_ptr<TMatrixRel> currentTMatRealm = getNewTMat();
  //  currentTMatRealm->evalNonRelMatrix(polePosEpsilonRealm);
  currentTMatRealm->evalMatrix(polePosEpsilonRealm, _orbitalL);

  std::shared_ptr<TMatrixRel> currentTMatImagp = getNewTMat();
  //  currentTMatImagp->evalNonRelMatrix(polePosEpsilonImagp);
  currentTMatImagp->evalMatrix(polePosEpsilonImagp, _orbitalL);

  std::shared_ptr<TMatrixRel> currentTMatImagm = getNewTMat();
  //  currentTMatImagm->evalNonRelMatrix(polePosEpsilonImagm);
  currentTMatImagm->evalMatrix(polePosEpsilonImagm, _orbitalL);

  InfoMsg << "\n\nm - i/2. Gamma: " << polePos.real()  << " - i/2. " << -2.*polePos.imag() << endmsg;

  for(unsigned int i=0 ; i<_phpVecs.size(); ++i) {
    ResidueProperties  currentResPropReal;
    ResidueProperties  currentResPropImag;
    ResidueProperties  currentResPropAverage;
    
    std::complex<double> resultEpsilonRealpInv= 1./(sqrt(_phpVecs.at(i)->factor(polePos, _orbitalL)) *
						    (*currentTMatRealp)(i,i)*sqrt(_phpVecs.at(i)->factor(polePos, _orbitalL)));
    std::complex<double> resultEpsilonRealmInv= 1./(sqrt(_phpVecs.at(i)->factor(polePos, _orbitalL)) *
						    (*currentTMatRealm)(i,i)*sqrt(_phpVecs.at(i)->factor(polePos, _orbitalL)));
    std::complex<double> resultApproxReal=(resultEpsilonRealpInv-resultEpsilonRealmInv)/(2.*epsilon);
    InfoMsg << "resultApproxRel: " << resultApproxReal << endmsg;
    
    std::complex<double> resultEpsilonImagpInv= 1./(sqrt(_phpVecs.at(i)->factor(polePos, _orbitalL)) *
						    (*currentTMatImagp)(i,i)*sqrt(_phpVecs.at(i)->factor(polePos, _orbitalL)));
    std::complex<double> resultEpsilonImagmInv= 1./(sqrt(_phpVecs.at(i)->factor(polePos, _orbitalL)) *
						    (*currentTMatImagm)(i,i)*sqrt(_phpVecs.at(i)->factor(polePos, _orbitalL)));
    std::complex<double> resultApproxImag=1./PawianConstants::i *
      (resultEpsilonImagpInv-resultEpsilonImagmInv)/(2.*epsilon);
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

    currentResPropReal.gammaigammaj=2.*abs(1./resultApproxReal)*2.*abs(1./resultApproxReal);
    currentResPropImag.gammaigammaj=2.*abs(1./resultApproxImag)*2.*abs(1./resultApproxImag);
    currentResPropAverage.gammaigammaj=2.*abs(1./resultApprox)*2.*abs(1./resultApprox);
    
    currentResPropReal.gammaiBRj=2.*abs(1./resultApproxReal)/(-2.*polePos.imag());
    currentResPropImag.gammaiBRj=2.*abs(1./resultApproxImag)/(-2.*polePos.imag());
    currentResPropAverage.gammaiBRj=2.*abs(1./resultApprox)/(-2.*polePos.imag());

      
    resPropReal.at(i)=currentResPropReal;
    resPropImag.at(i)=currentResPropImag;
    resPropAverage.at(i)=currentResPropAverage;
  }
}

void TMatrixResidueExtr::printResults(std::complex<double> polePos, std::vector<ResidueProperties> resPropAv){
  InfoMsg << "******* Summary Residues *******" << endmsg;
  InfoMsg << "M - i/2 Gamma = " << polePos.real()  << " - i/2 " << -2.*polePos.imag() << endmsg;

  for (unsigned int i=0; i<resPropAv.size(); ++i){
    InfoMsg << i
	    << "\t" << resPropAv.at(i).gammai << " +/- " << resPropAv.at(i).errGammai
	    << "\t" << -resPropAv.at(i).gammai/(2.*polePos.imag()) <<  " +/- " << -resPropAv.at(i).errGammai/(2.*polePos.imag())
	    << endmsg;
  }
}

void TMatrixResidueExtr::dumpResult(std::complex<double> polePos, std::vector<ResidueProperties> resPropReal, 
				    std::vector<ResidueProperties> resPropImag, 
				    std::vector<ResidueProperties> resPropAv){
  std::string oFileName="residues_"+_extractionMethod+".out";
  std::ofstream theStream(oFileName);
  double gammaTotalr=0.;
  double gammaTotali=0.;
  double gammaTotalav=0.;
  double gammaTotalavError=0.;
  double BRTotalr=0.;
  double BRTotali=0.;
  double BRTotalav=0.;
  double BRTotalavError=0.;

  theStream << "m - i/2. Gamma = " << polePos.real()  << " - i/2. " << -2.*polePos.imag() << std::endl << std::endl;
  for (unsigned int i=0; i<resPropReal.size(); ++i){
    theStream << "\nchannel :              " << i 
              << "\nTheta_r :              "  << resPropReal.at(i).theta*180./M_PI 
	      << "\nTheta_i :              "  << resPropImag.at(i).theta*180./M_PI 
	      << "\nTheta_av:              "  << resPropAv.at(i).theta*180./M_PI 
              << "\nGamma_i(r) : " << resPropReal.at(i).gammai
              << "\nGamma_i(i) : " << resPropImag.at(i).gammai
      	      << "\nGamma_i(av): " << resPropAv.at(i).gammai
	      << " +/- " << resPropAv.at(i).errGammai
	      << "\nBR[%]_i(r)   : " << resPropReal.at(i).gammaiBRj 
	      << "\nBR[%]_i(i)   : " << resPropImag.at(i).gammaiBRj 
	      << "\nBR[%]_i(av)   : " << resPropAv.at(i).gammaiBRj
	      <<  " +/- " << resPropAv.at(i).errGammaiBRj
              << std::endl;

        gammaTotalr+=resPropReal.at(i).gammai;
	gammaTotali+=resPropImag.at(i).gammai;
	gammaTotalav+=resPropAv.at(i).gammai;
	gammaTotalavError+=resPropAv.at(i).errGammai*resPropAv.at(i).errGammai;
	BRTotalr+=resPropReal.at(i).gammaiBRj;
	BRTotali+=resPropImag.at(i).gammaiBRj;
	BRTotalav+=resPropAv.at(i).gammaiBRj;
	BRTotalavError+=resPropAv.at(i).errGammaiBRj;	
  }

  gammaTotalavError=std::sqrt(gammaTotalavError);
  BRTotalavError=std::sqrt(BRTotalavError);  
  theStream << "\nall" 
	    << "\nGamma_total(r): " << gammaTotalr
	    << "\nGamma_total(i): " << gammaTotali
	    << "\nGamma_total   : " << gammaTotalav
	    << " +/- " << gammaTotalavError
	    << "\n\nBR_total(r): " << BRTotalr
	    << "\nBR_total(i): " << BRTotali
	    << "\nBR_total   : " << BRTotalav
	    << " +/- " << BRTotalavError    
	    << std::endl;
}

