//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
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
//#include <boost/multi_array.hpp>
#include "Examples/Tutorial/LineShapes/TMatrixGeneral.hh"
#include "Examples/Tutorial/LineShapes/RiemannSheetAnalyzer.hh"
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
#include "Utils/PawianConstants.hh"

#include "ConfigParser/ParserBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/TMatrixDynamics.hh"
#include "FitParams/ParamFactory.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

TMatrixGeneral::TMatrixGeneral(std::string pathToConfigParser, std::string pathToFitParams, int numStepsForSheetScan, std::vector<double> energyPlaneBorders) :
  _noOfSteps(500)
  ,_stepSize(0.)
  ,_massMin(100000.)
  ,_massMax(0.)
  ,_kMatrixParser(new KMatrixParser(pathToConfigParser))
  ,_theTFile(0)
  ,_pathToFitParams(pathToFitParams)
  ,_orbitalL(0)
{ 
  init();
  std::string rootFileName="./TMatrixGeneral.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");

  _massMax+=0.4;
  _stepSize=(_massMax-_massMin)/_noOfSteps;

  std::string ampRealKey="AmpReal";
  std::string ampImagKey="AmpImag";
  std::string argandKey="Argand";
  std::string phaseKey="Phase";
  std::string inelasticityKey="Elasticity";
  std::string ImagT11m1Key="Imag(T_{11}^{-1}";

  for(unsigned int i=0; i<_gFactorNames.size(); ++i){
    std::string key=_gFactorNames.at(i);
    std::string currentAmpRealKey=ampRealKey+" "+key;
    TH1F* currentAmpRealH1=new TH1F(currentAmpRealKey.c_str(),currentAmpRealKey.c_str(), _noOfSteps-1, _massMin, _massMax);
    currentAmpRealH1->SetYTitle(currentAmpRealKey.c_str());
    currentAmpRealH1->SetXTitle("mass");
    _AmpRealH1Vec.push_back(currentAmpRealH1);

    std::string currentAmpImagKey=ampImagKey+" "+key;
    TH1F* currentAmpImagH1=new TH1F(currentAmpImagKey.c_str(),currentAmpImagKey.c_str(), _noOfSteps-1, _massMin, _massMax);
    currentAmpImagH1->SetYTitle(currentAmpImagKey.c_str());
    currentAmpImagH1->SetXTitle("mass");
    _AmpImagH1Vec.push_back(currentAmpImagH1);


    std::string currentImagT11m1Key=ImagT11m1Key+" "+key;
    TH1F* currentImagT11m1H1=new TH1F(currentImagT11m1Key.c_str(), currentImagT11m1Key.c_str(), _noOfSteps-1, _massMin, _massMax);
    currentImagT11m1H1->SetYTitle(currentImagT11m1Key.c_str());
    currentImagT11m1H1->SetXTitle("mass");
    _ImagT11m1H1Vec.push_back(currentImagT11m1H1);

    std::string currentArgandKey=argandKey+key;
    TH2F* currentArgandH2=new TH2F(currentArgandKey.c_str(), currentArgandKey.c_str(), 200, -1., 1., 200, -.5, 1.5);
    std::string currentArgandYKey="Imag"+key;
    std::string currentArgandXKey="Real"+key;
    currentArgandH2->SetYTitle(currentArgandYKey.c_str());
    currentArgandH2->SetXTitle(currentArgandXKey.c_str());
    _ArgandH2Vec.push_back(currentArgandH2);

    std::string currentPhaseKey=phaseKey+key;
    TH2F* currentPhaseH2=new TH2F(currentPhaseKey.c_str(), currentPhaseKey.c_str(), _noOfSteps-1, _massMin, _massMax, 360, 0., 360.);

    std::string currentPhaseYKey="#delta"+key;
    currentPhaseH2->SetYTitle(currentPhaseYKey.c_str());
    currentPhaseH2->SetXTitle("mass");
    _PhaseH2Vec.push_back(currentPhaseH2);

    std::string currentElasticityKey=inelasticityKey+key;
    TH1F* currentElasticityH1=new TH1F(currentElasticityKey.c_str(), currentElasticityKey.c_str(), _noOfSteps-1, _massMin, _massMax);
    std::string currentElasticityYKey="#eta"+key;
    currentElasticityH1->SetYTitle(currentElasticityYKey.c_str());
    currentElasticityH1->SetXTitle("mass/GeV");
    _ElasticityH1Vec.push_back(currentElasticityH1);

    std::string currentT11Key="|T_{11}|^{2}"+key;
    TH1F* currentSqrT11H1=new TH1F(currentT11Key.c_str(), currentT11Key.c_str(), _noOfSteps-1, _massMin, _massMax);
    currentSqrT11H1->SetYTitle("|T_{11}|^{2}");
    currentSqrT11H1->SetXTitle("mass/GeV");
    _SqrT11H1Vec.push_back(currentSqrT11H1);

    std::string currentphpKey="phase space factor"+key;
    TH1F* currentphpH1=new TH1F(currentphpKey.c_str(), currentphpKey.c_str(), _noOfSteps-1, _massMin, _massMax);
    currentphpH1->SetYTitle("#rho");
    currentphpH1->SetXTitle("mass/GeV");
    _phpH1Vec.push_back(currentphpH1);

    std::string currentphpRealKey="phase space factor (real)"+key;
    TH1F* currentphpRealH1=new TH1F(currentphpRealKey.c_str(), currentphpRealKey.c_str(), _noOfSteps-1, _massMin, _massMax);
    currentphpRealH1->SetYTitle("Re #rho");
    currentphpRealH1->SetXTitle("mass/GeV");
    _phpH1RealVec.push_back(currentphpRealH1);

    std::string currentphpImagKey="phase space factor (imag)"+key;
    TH1F* currentphpImagH1=new TH1F(currentphpImagKey.c_str(), currentphpImagKey.c_str(), _noOfSteps-1, _massMin, _massMax);
    currentphpImagH1->SetYTitle("Imag #rho");
    currentphpImagH1->SetXTitle("mass/GeV");
    _phpH1ImagVec.push_back(currentphpImagH1);

    std::string currentSqrtT1iKey="sqrtT1i"+key;
    TH1F* currentSqrtT1iH1=new TH1F(currentSqrtT1iKey.c_str(), currentSqrtT1iKey.c_str(), _noOfSteps-1, _massMin, _massMax);
    currentSqrtT1iH1->SetYTitle("#sqrt{#rho_{1}}#sqrt{#rho_{i}} |T_{1i}|");
    currentSqrtT1iH1->SetXTitle("mass/GeV");
    _SqrtT1iH1Vec.push_back(currentSqrtT1iH1);

    std::string currentArgandUnits1iKey="ArgandUnit1i"+key;
    TH1F* currentArgandUnits1iH1=new TH1F(currentArgandUnits1iKey.c_str(), currentArgandUnits1iKey.c_str(), _noOfSteps-1, _massMin, _massMax);
    currentArgandUnits1iH1->SetYTitle("|#sqrt{#rho_{1}}#sqrt{#rho_{i}} T_{1i}|");
    currentArgandUnits1iH1->SetXTitle("mass/GeV");
    _ArgandUnits1iH1Vec.push_back(currentArgandUnits1iH1);

    std::string currentSqrS1iKey="sqrS1i"+key;
    TH1F* currentSqrS1iH1=new TH1F(currentSqrS1iKey.c_str(), currentSqrS1iKey.c_str(), _noOfSteps-1, _massMin, _massMax);
    currentSqrS1iH1->SetYTitle("|#sqrt{#rho_{1}}#sqrt{#rho_{i}} S_{1i}|^{2}");
    currentSqrS1iH1->SetXTitle("mass/GeV");
    _SqrS1iH1Vec.push_back(currentSqrS1iH1);

    std::string currentdelta1iKey="delta1i"+key;
    TH1F* currentDelta1iH1=new TH1F(currentdelta1iKey.c_str(), currentdelta1iKey.c_str(), _noOfSteps-1, _massMin, _massMax);
    currentDelta1iH1->SetYTitle("#delta1i/grad");
    currentDelta1iH1->SetXTitle("mass/GeV");
    _delta1iVec.push_back(currentDelta1iH1);

    std::string currentSppediiKey="speedii"+key;
    TH1F* currentSpeediiH1=new TH1F(currentSppediiKey.c_str(), currentSppediiKey.c_str(), _noOfSteps-1, _massMin, _massMax);
    currentSpeediiH1->SetYTitle("speed");
    currentSpeediiH1->SetXTitle("mass/GeV");
    _speedPlotH1Vec.push_back(currentSpeediiH1);
  }

  DebugMsg << "_massMin: "<< _massMin
  	   << "\n_massMax: "<< _massMax
  	   << "\n_stepSize: " << _stepSize
  	   << endmsg;

  std::ofstream oStream;
  oStream.open("scatteringOut.txt");
  std::ofstream oStream_phi;
  oStream_phi.open("scatteringOut_phi.txt");
 double oldT00RelReal=1.;
  int n180ShiftRel(0);

  double massOld=0.;
  vector< complex<double> > TiiOld;
  TiiOld.resize(_gFactorNames.size());
 
  for (double mass=_massMin+_stepSize/0.5; mass<_massMax; mass+=_stepSize){
    Vector4<double> mass4Vec(mass, 0.,0.,0.);
    _tMatr->evalMatrix(mass);

   std::shared_ptr<TMatrixRel> tMatrInv=std::shared_ptr<TMatrixRel>(new TMatrixRel(_kMatr));
   tMatrInv->evalMatrix(mass);
   tMatrInv->invert();
   double delta1=0.;

    for(unsigned int i=0; i<_gFactorNames.size(); ++i){
      complex<double> currentRho=_phpVecs.at(i)->factor(mass);

      //      if(currentRho.real() > 1.e-8){ 
	// _AmpRealH1Vec.at(i)->Fill(mass, sqrt(currentRho.real())*(*_tMatr)(i,i).real());
	// _AmpImagH1Vec.at(i)->Fill(mass, sqrt(currentRho.real())*(*_tMatr)(i,i).imag());
	_AmpRealH1Vec.at(i)->Fill(mass, (*_tMatr)(i,i).real());
	_AmpImagH1Vec.at(i)->Fill(mass, (*_tMatr)(i,i).imag());	
	_ImagT11m1H1Vec.at(i)->Fill(mass, -(*tMatrInv)(i,i).imag());
	
	_ArgandH2Vec.at(i)->Fill( currentRho.real()*(*_tMatr)(i,i).real(), currentRho.real()*(*_tMatr)(i,i).imag());
	double currentphase=180.*atan2((*_tMatr)(i,i).imag(),(*_tMatr)(i,i).real()) / 3.1415;
	complex<double> SijRel_phi=complex<double>(1.,0.)+2.*PawianConstants::i*currentRho.real()*(*_tMatr)(i,i);
  if(i==0) {
	oStream_phi << mass << "\t" 
			<< currentphase << "\t" << 0.01 << "\t" 
		<<  sqrt(norm(SijRel_phi)) << "\t" << 0.001 
			<< endl;
  }
	_PhaseH2Vec.at(i)->Fill(mass, currentphase);
	// double sqrtFactor=(*_tMatr)(i,i).real()*(*_tMatr)(i,i).real()+((*_tMatr)(i,i).imag()-0.5)*((*_tMatr)(i,i).imag()-0.5);
	// double currentElasticity=2.*sqrt(sqrtFactor);
	//      }
      complex<double> S00_rel=complex<double>(1.,0.)+2.*complex<double>(0.,1.)*currentRho.real()*(*_tMatr)(i,i);
      _ElasticityH1Vec.at(i)->Fill(mass, sqrt(norm(S00_rel)));

      _SqrT11H1Vec.at(i)->Fill(mass,currentRho.real()*norm((*_tMatr)(i,i)));

      _phpH1Vec.at(i)->Fill(mass, sqrt(norm(currentRho)));
      _phpH1RealVec.at(i)->Fill(mass, currentRho.real());
      _phpH1ImagVec.at(i)->Fill(mass, currentRho.imag());

      _SqrtT1iH1Vec.at(i)->Fill(mass, sqrt(_phpVecs.at(0)->factor(mass).real()*currentRho.real())*sqrt(norm((*_tMatr)(0,i))) );

      _ArgandUnits1iH1Vec.at(i)->Fill(mass, (2.*_orbitalL+1.)*norm(sqrt(_phpVecs.at(0)->factor(mass).real()*currentRho.real())*(*_tMatr)(0,i)));

      complex<double> currentTii=currentRho.real()*(*_tMatr)(i,i);

    if (mass > _massMin+_stepSize){
      complex<double> deltaTii = currentTii-TiiOld.at(i);
      double deltaMass=mass-massOld;
      double speed = std::abs(deltaTii/deltaMass);
      //    InfoMsg << "speed(" << mass << "): " << speed << endmsg;
      // InfoMsg << "deltaT1i: " << deltaT1i << endmsg;
      // InfoMsg << "deltaMass: " << deltaMass << endmsg;
      _speedPlotH1Vec.at(i)->Fill(mass,speed); 
    }
    
    TiiOld[i]=currentTii;

      if(currentRho.real()>0.){
	complex<double> S0iRel=2.*PawianConstants::i*sqrt(_phpVecs.at(0)->factor(mass).real()*currentRho.real())*(*_tMatr)(0,i);
	if(i==0) S0iRel+=complex<double>(1.,0.);  
	_SqrS1iH1Vec.at(i)->Fill(mass, norm(S0iRel));
      }

      complex<double> currentTiiRel_rho= (*_tMatr)(i,i)*currentRho.real();
      //complex<double> currentTiiRel_rho= (*_tMatr)(i,i);
      double currentReEiiRel = currentTiiRel_rho.real();
      double currentImEiiRel = currentTiiRel_rho.imag() - 0.5;
      // InfoMsg << "currentReEiiRel: " << currentReEiiRel << "\tcurrentImEiiRel: " << currentImEiiRel << endmsg;
      double deltaiiRel = 0.5*atan2(currentImEiiRel, fabs(currentReEiiRel))*PawianConstants::radToDeg + 45.0;
      if (currentReEiiRel  < 0.0) {deltaiiRel = 180.0 - deltaiiRel;}
      if(i==0){
	delta1=deltaiiRel;
	_delta1iVec.at(i)->Fill(mass, delta1);
      }
      else if(currentRho.real() > 1.e-8) _delta1iVec.at(i)->Fill(mass, delta1+deltaiiRel);

      if (i==0){
	complex<double> currentTijRel_rho= (*_tMatr)(i,i)*currentRho.real();
	double currentReE = currentTijRel_rho.real();
	double currentImE = currentTijRel_rho.imag()- 0.5;
	double theDelta = 0.5*atan2( currentImE, fabs(currentReE) )*PawianConstants::radToDeg + 45.0;
	if (currentReE   < 0.0) {theDelta = 180.0 - theDelta;}
	//double theDelta = 0.5*atan2( currentImE, currentReE )*PawianConstants::radToDeg;
	 // Have we gone through 180 deg (or 2*delta through 360 deg)?
	 // if (oldT00RelReal < 0.0 && currentReE > 0.0) {n180ShiftRel = 1;}
	theDelta += 180.0*n180ShiftRel;
	if (oldT00RelReal < 0.0 && currentReE > 0.0) {n180ShiftRel += 1;}
	// if ( oldDelta > 1.e-8 && fabs(theDelta-oldDelta) > 90.) {
	//   n180ShiftRel += 1;
	// }
	//	 theDelta += 180.0*n180ShiftRel; 
	  //	if (oldT00RelReal < 0.0 && currentReE > 0.0) {n180ShiftRel += 1;}
	// thePhitheDelta += 180.0*n180ShiftRel;
	//theDelta += 180.0*n180ShiftRel;
	//	double thePhi=0.5*currentphase+ 45.0;
	complex<double> SijRel=complex<double>(1.,0.)+2.*PawianConstants::i*currentRho.real()*(*_tMatr)(i,i);
  if(abs((*_tMatr)(i,i))>0.1) {
	// theDelta=360.*atan2((*_tMatr)(i,i).imag(),(*_tMatr)(i,i).real()) / 3.1415;
	oStream << mass << "\t" 
			<< theDelta << "\t" << 0.01 << "\t" 
		<<  sqrt(norm(SijRel)) << "\t" << 0.001 
			<< endl;
  }
	oldT00RelReal = currentReE;
	//        oldDelta=theDelta; 
      }
    }
    massOld=mass;
  }
  
  oStream.close();
  oStream_phi.close();
  
  if(energyPlaneBorders[0] == 0)
    energyPlaneBorders[0] = _massMin;
  if(energyPlaneBorders[2] == 0)
    energyPlaneBorders[2] = _massMax;
  if(energyPlaneBorders[1] > 0.) energyPlaneBorders[1]=-energyPlaneBorders[1]; //minImagMass should be nagative
  //search projection index

  const std::string porjectionParticleNames=_kMatrixParser->projection();
  std::istringstream projParticles(porjectionParticleNames);
  std::string firstProjParticleName;
  std::string secondProjParticleName;
  projParticles >> firstProjParticleName >> secondProjParticleName;
  std::string projKey=firstProjParticleName+secondProjParticleName;

  unsigned int projectionIndex=0;   
  bool found=false;
  for(unsigned int i=0; i<_gFactorNames.size();++i){
    if(projKey==_gFactorNames[i]){
      projectionIndex=i;
      found=true;
    }
  }

  if (!found){
    Alert << "projection index for key " << projKey << " not found" << endmsg;
    exit(0);
  }

  // RiemannSheetAnalyzer(_kMatrixParser->noOfChannels(), _tMatr,
  // 		       std::complex<double>(energyPlaneBorders[0], energyPlaneBorders[1]),
  // 		       std::complex<double>(energyPlaneBorders[2], energyPlaneBorders[3]),
  // 		       numStepsForSheetScan, 
  // 		       projectionIndex);
}

TMatrixGeneral::~TMatrixGeneral()
{
  _theTFile->Write();
  _theTFile->Close();
}


void TMatrixGeneral::init(){
  GlobalEnv::instance()->setup();
  _particleTable=GlobalEnv::instance()->particleTable();
  std::shared_ptr<TMatrixDynamics> tMatDynPtr= std::shared_ptr<TMatrixDynamics>(new TMatrixDynamics(_kMatrixParser));

    std::shared_ptr<AbsPawianParameters> params=ParamFactory::instance()->getParametersPointer("Pawian");
    tMatDynPtr->fillDefaultParams(params);
 
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

  InfoMsg << "The k-Matrix input parameter are: " << endmsg;
  _params->print(std::cout);
  if(_pathToFitParams != "") tMatDynPtr->updateFitParams(_params); 

  _kMatr = tMatDynPtr->getKMatix(); 
  _tMatr = tMatDynPtr->getTMatix();  

  _phpVecs=_kMatr->phaseSpaceVec();
  _gFactorNames= tMatDynPtr->gFactorNames();
  _orbitalL= tMatDynPtr->orbitalL();

  std::vector<std::string> poleNameAndMassVecs=_kMatrixParser->poles();
  std::vector<std::string>::iterator itString;
  for (itString=poleNameAndMassVecs.begin(); itString!=poleNameAndMassVecs.end(); ++itString){
    std::istringstream poleIString(*itString);
    std::string currentPoleName;
    std::string currentPoleMassStr;
    poleIString >> currentPoleName >> currentPoleMassStr;
    
    std::istringstream currentPoleMassiStr(currentPoleMassStr);
    double currentValue;
    if(!(currentPoleMassiStr >> currentValue)){
      Alert << "cannot convert " << currentPoleMassStr << " to a double value" << endmsg;
      exit(0);
    }
    if (currentValue>_massMax) _massMax=currentValue;
    
  }

  const std::vector<std::string> gFacStringVec=_kMatrixParser->gFactors();
  DebugMsg << "gFacStringVec.size(): " << gFacStringVec.size() << endmsg;
  std::map<std::pair<std::string, std::string>, std::string> phpDescriptionVec=_kMatrixParser->phpDescriptionMap();
  std::cout << "phpDescriptionVec.size(): " << phpDescriptionVec.size() << std::endl;
  std::vector<std::string>::const_iterator itStrConst;
  for(itStrConst=gFacStringVec.begin(); itStrConst!=gFacStringVec.end(); ++itStrConst){
    std::istringstream particles(*itStrConst);
    std::string firstParticleName;
    std::string secondParticleName;
    particles >> firstParticleName >> secondParticleName;
    std::pair<std::string, std::string> currentParticlePair=make_pair(firstParticleName, secondParticleName);
    Particle* firstParticle = _particleTable->particle(firstParticleName);
    Particle* secondParticle = _particleTable->particle(secondParticleName);
    if(0==firstParticle || 0==secondParticle){
      Alert << "particle with name: " << firstParticleName <<" or " << secondParticleName << " doesn't exist in pdg-table" << endmsg;
      exit(0);
    }
    double currentMassSum=firstParticle->mass()+secondParticle->mass();
    if(currentMassSum<_massMin) _massMin=currentMassSum;
  }

}

