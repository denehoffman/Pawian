//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)  
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
  _pathToSerialzationFile(pathToSerialzationFile)
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

void TMatrixResidueExtr::GetCovMatrix(){
  std::ifstream serializationStream(_pathToSerialzationFile.c_str());

  if(!serializationStream.is_open()){
	Alert << "Could not open serialization file." << endmsg;
	exit(0);
  }
  _thePwaCovMatrix = std::shared_ptr<PwaCovMatrix>(new PwaCovMatrix);
  boost::archive::text_iarchive boostInputArchive(serializationStream);

  boostInputArchive >> *_thePwaCovMatrix;
}

std::complex<double> TMatrixResidueExtr::CalcMassWidth(){
  //  _tMatFit->updateTMatDy(theFitParams);
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

void  TMatrixResidueExtr::CalcResidueAll(){
  std::complex<double> result(0.,0.);
  const double epsilon=0.000001;
  std::complex<double> polePos = CalcMassWidth();
  std::complex<double> polePosEpsilonImagp = polePos + std::complex<double>(0., epsilon);
  std::complex<double> polePosEpsilonImagm = polePos + std::complex<double>(0., -epsilon);
  std::complex<double> polePosEpsilonRealp = polePos + std::complex<double>(epsilon, 0.);
  std::complex<double> polePosEpsilonRealm = polePos + std::complex<double>(-epsilon, 0.);
  // std::complex<double> polePosEpsilonImagp(-polePos.imag()+epsilon, -polePos.real()+epsilon);
  // std::complex<double> polePosEpsilonImagm(-polePos.imag()-epsilon, -polePos.real()-epsilon);
  // std::complex<double> polePosEpsilonRealp(polePos.real()+epsilon, polePos.imag()+epsilon);
  // std::complex<double> polePosEpsilonRealm(polePos.real()-epsilon, polePos.imag()-epsilon);
  
  std::shared_ptr<TMatrixRel> currentTMatRealp=_tMatFit->getNewTMat();
  currentTMatRealp->evalMatrix(polePosEpsilonRealp);
  // currentTMatRealp->invert();

  std::shared_ptr<TMatrixRel> currentTMatRealm=_tMatFit->getNewTMat();
  currentTMatRealm->evalMatrix(polePosEpsilonRealm);
  // currentTMatRealm->invert();

  std::shared_ptr<TMatrixRel> currentTMatImagp=_tMatFit->getNewTMat();
  currentTMatImagp->evalMatrix(polePosEpsilonImagp);
  // currentTMatImagp->invert();

  std::shared_ptr<TMatrixRel> currentTMatImagm=_tMatFit->getNewTMat();
  currentTMatImagm->evalMatrix(polePosEpsilonImagm);
  // currentTMatImagm->invert();

  InfoMsg << "\n\nm - i/2. Gamma: " << polePos.real()  << " - i/2. " << -2.*polePos.imag() << endmsg;
  double gammaSum=0.;
  double brSum=0.;  
  for(int i=0 ; i<_phpVecs.size(); ++i){
    std::complex<double> resultEpsilonRealpInv= 1./(*currentTMatRealp)(i,i);
    std::complex<double> resultEpsilonRealmInv= 1./(*currentTMatRealm)(i,i);
    // std::complex<double> resultEpsilonRealpInv= (*currentTMatRealp)(i,i);
    // std::complex<double> resultEpsilonRealmInv= (*currentTMatRealm)(i,i);        
    std::complex<double> resultApproxReal=(resultEpsilonRealpInv-resultEpsilonRealmInv)/(2.*epsilon);
    InfoMsg << "resultApproxRel: " << resultApproxReal << endmsg;

    std::complex<double> resultEpsilonImagpInv= 1./(*currentTMatImagp)(i,i);
    std::complex<double> resultEpsilonImagmInv= 1./(*currentTMatImagm)(i,i);
    // std::complex<double> resultEpsilonImagpInv= (*currentTMatImagp)(i,i);
    // std::complex<double> resultEpsilonImagmInv= (*currentTMatImagm)(i,i);
    std::complex<double> resultApproxImag=1./PawianConstants::i*(resultEpsilonImagpInv-resultEpsilonImagmInv)/(2.*epsilon);
    InfoMsg << "resultApproxImag: " << resultApproxImag << endmsg;

    std::complex<double> resultApprox = (resultApproxReal+resultApproxImag)/2.;
    //    std::complex<double> resultApprox = resultApproxReal;
    InfoMsg << "resultApprox: " << resultApprox << endmsg;
    
    double rhoi=imag(_phpVecs.at(i)->ChewM(polePos));
    double resAbsRi=abs(1./resultApprox*rhoi);
    double resThetai=-M_PI+atan(imag(1./resultApprox*rhoi)/real(1./resultApprox*rhoi));

    InfoMsg << "\nrho(" << i << "): " << rhoi << endmsg; 
    InfoMsg << "\nresAbsR(" << i << "): " << resAbsRi << endmsg; 
    InfoMsg << "resTheta(" << i << "): " << resThetai << " rad\t" << resThetai*180./M_PI << " deg" << endmsg;

    double g_i=sqrt(resAbsRi);
    //double gammai=calcPartialWidth(g_i, polePos, _phpVecs.at(i));
    double gammai=calcPartialWidth(1./resultApprox, polePos, _phpVecs.at(i));
    //if(gammai>0.) gammaSum+=gammai;
    gammaSum+=fabs(gammai);
    InfoMsg << "Gamma(" << i << "): " << gammai << endmsg;
    InfoMsg << "BR(" << i << "): " << gammai/(-2.*polePos.imag()) << endmsg;
    //    if(gammai>0.) brSum+=gammai/(-2.*polePos.imag());
    brSum+=fabs(gammai/(-2.*polePos.imag()));     
  }
  InfoMsg << "GammaSum: " << gammaSum << "\tBR Sum: " << brSum << endmsg;
}

double TMatrixResidueExtr::calcPartialWidth(double gFac, std::complex<double> poleMass, std::shared_ptr<AbsPhaseSpace> php){
  double result=0.;
  std::shared_ptr<PhaseSpaceIsobar> phpIso =  std::dynamic_pointer_cast<PhaseSpaceIsobar>(php);

  double decMass1=phpIso->mass1();
  double decMass2=phpIso->mass2();

  complex<double> breakUmMom = php->breakUpMom(poleMass);
  //  complex<double> breakUmMom = php->breakUpMom(poleMass.real());
  //  double breakUmMom = real(php->breakUpMom(poleMass));  
  double energy1 = sqrt(breakUmMom.real()*breakUmMom.real()+decMass1*decMass1);
  double energy2 = sqrt(breakUmMom.real()*breakUmMom.real()+decMass2*decMass2);

  double gTildeFac = 2.*M_PI*gFac*sqrt(energy1*energy2/decMass1);
  result=gTildeFac*gTildeFac/(2.*M_PI)*decMass1*breakUmMom.real()/poleMass.real();
  result*=2.; //result seems to be Gamma/2
  // double decMass12=php->thresholdMass()/2.;
  // complex<double> breakUmMom = php->breakUpMom(poleMass);
  // //  double breakUmMom = real(php->breakUpMom(poleMass));  
  // double energy12 = sqrt(breakUmMom.real()*breakUmMom.real()+decMass12*decMass12);
  // double gTildeFac = 2.*M_PI*gFac*sqrt(energy12*energy12/decMass12);
  // result=gTildeFac*gTildeFac/(2.*M_PI)*decMass12*breakUmMom.real()/poleMass.real();
  
  return result;  
}


double  TMatrixResidueExtr::calcPartialWidth(std::complex<double> am1, std::complex<double> poleMass, std::shared_ptr<AbsPhaseSpace> php){
  double result=0.;
  //  complex<double> am1Mod(am1.real(), 2.*am1.imag());
  complex<double> gFac=sqrt(am1);
  //  complex<double> gFac=sqrt(am1);
  std::shared_ptr<PhaseSpaceIsobar> phpIso =  std::dynamic_pointer_cast<PhaseSpaceIsobar>(php);

  double decMass1=phpIso->mass1();
  double decMass2=phpIso->mass2();
  double breakUmMom = real(php->breakUpMom(poleMass));
  double energy1 = sqrt(breakUmMom*breakUmMom+decMass1*decMass1);
  double energy2 = sqrt(breakUmMom*breakUmMom+decMass2*decMass2);

  complex<double> gTildeFac = 2.*M_PI*gFac*sqrt(energy1*energy2/decMass1);

  result=norm(gTildeFac)/(2.*M_PI)*decMass1*breakUmMom/poleMass.real();
  //  result=norm(gTildeFac)/(2.*M_PI)*decMass1/poleMass.real();
    
  InfoMsg << "\ndecMass1: " <<  decMass1 << "\tdecMass2: " <<  decMass2 << endmsg;
  InfoMsg << "breakUmMom: " <<  breakUmMom << endmsg;
  InfoMsg << "energy1: " <<  energy1 << "\tenergy2: " << energy2  << endmsg;

  result/=0.5*M_PI*energy1*energy2;
  //  result*=0.2*M_PI/energy2;
  //  result*=decMass1*M_PI/energy2;  
  return result;
}
