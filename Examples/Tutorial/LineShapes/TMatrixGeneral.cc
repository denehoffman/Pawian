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

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

TMatrixGeneral::TMatrixGeneral(std::string pathToConfigParser, int numStepsForSheetScan, std::vector<double> energyPlaneBorders) :
  _noOfSteps(1000)
  ,_stepSize(0.)
  ,_massMin(100000.)
  ,_massMax(0.)
  ,_kMatrixParser(new KMatrixParser(pathToConfigParser))
  ,_theTFile(0)
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
  }

  DebugMsg << "_massMin: "<< _massMin
  	   << "\n_massMax: "<< _massMax
  	   << "\n_stepSize: " << _stepSize
  	   << endmsg;

  for (double mass=_massMin+_stepSize/0.5; mass<_massMax; mass+=_stepSize){
    Vector4<double> mass4Vec(mass, 0.,0.,0.);
    _tMatr->evalMatrix(mass);

   std::shared_ptr<TMatrixRel> tMatrInv=std::shared_ptr<TMatrixRel>(new TMatrixRel(_kMatr));
   tMatrInv->evalMatrix(mass);
   tMatrInv->invert();

    for(unsigned int i=0; i<_gFactorNames.size(); ++i){
      complex<double> currentRho=_phpVecs.at(i)->factor(mass);

      _AmpRealH1Vec.at(i)->Fill(mass, sqrt(currentRho.real())*(*_tMatr)(i,i).real());
      _AmpImagH1Vec.at(i)->Fill(mass, sqrt(currentRho.real())*(*_tMatr)(i,i).imag());
      _ImagT11m1H1Vec.at(i)->Fill(mass, -(*tMatrInv)(i,i).imag());
       
      _ArgandH2Vec.at(i)->Fill( currentRho.real()*(*_tMatr)(i,i).real(), currentRho.real()*(*_tMatr)(i,i).imag());
      double currentphase=360.*atan2((*_tMatr)(i,i).imag(),(*_tMatr)(i,i).real()) / 3.1415;
      _PhaseH2Vec.at(i)->Fill(mass, currentphase);
      // double sqrtFactor=(*_tMatr)(i,i).real()*(*_tMatr)(i,i).real()+((*_tMatr)(i,i).imag()-0.5)*((*_tMatr)(i,i).imag()-0.5);
      // double currentElasticity=2.*sqrt(sqrtFactor);
      complex<double> S00_rel=complex<double>(1.,0.)+2.*complex<double>(0.,1.)*currentRho.real()*(*_tMatr)(i,i);
      _ElasticityH1Vec.at(i)->Fill(mass, sqrt(norm(S00_rel)));

      _SqrT11H1Vec.at(i)->Fill(mass,currentRho.real()*norm((*_tMatr)(i,i)));

      _phpH1Vec.at(i)->Fill(mass, sqrt(norm(currentRho))); 
    }    
  }

  if(energyPlaneBorders[0] == 0)
    energyPlaneBorders[0] = _massMin;
  if(energyPlaneBorders[2] == 0)
    energyPlaneBorders[2] = _massMax;

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

  RiemannSheetAnalyzer(_kMatrixParser->noOfChannels(), _tMatr,
		       std::complex<double>(energyPlaneBorders[0], energyPlaneBorders[1]),
		       std::complex<double>(energyPlaneBorders[2], energyPlaneBorders[3]),
		       numStepsForSheetScan, 
		       projectionIndex);
}

TMatrixGeneral::~TMatrixGeneral()
{
  _theTFile->Write();
  _theTFile->Close();
}


void TMatrixGeneral::init(){
   PdtParser pdtParser;
   std::string theSourcePath=getenv("TOP_DIR");
   std::string pdtFileRelPath="/Particle/pdtNew.table";
   std::string pdtFile(theSourcePath+pdtFileRelPath);
   _particleTable = new ParticleTable;

   if (!pdtParser.parse(pdtFile, *_particleTable)) {
      Alert << "can not parse particle table " << pdtFile << endmsg;
      exit(1);
   }


  std::vector< std::string> poleNames;
  std::vector<double> currentPoleMasses;

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
    
    _poleMasses.push_back(currentValue);
    _poleNames.push_back(currentPoleName);
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

    std::string currentPhpDescription = phpDescriptionVec.at(currentParticlePair);

    std::vector<double> fsParticleMasses;
    fsParticleMasses.push_back(firstParticle->mass());
    fsParticleMasses.push_back(secondParticle->mass()); 
   
    std::shared_ptr<AbsPhaseSpace> currentPhp=PhaseSpaceFactory::instance()->getPhpPointer(currentPhpDescription, fsParticleMasses);

    _phpVecs.push_back(currentPhp);

    std::string gFactorKey=firstParticleName+secondParticleName;    
    _gFactorNames.push_back(gFactorKey); 
    

    for (int i=0; i<int(_kMatrixParser->noOfPoles()); ++i){
      std::string currentPoleName=_poleNames[i];
      std::string currentgValueStr;
      if(!(particles >> currentgValueStr)){
        Alert << "g-factors for pole " << currentPoleName << " does not exist!" << endmsg;
        exit(0);
      }
      std::istringstream currentgValueiStr(currentgValueStr);
      double currentGValue;
      if (!(currentgValueiStr >> currentGValue)){
        Alert << "cannot convert " << currentgValueStr << " to a double value" << endmsg;
        exit(0);
      }
      std::string gFactorKey=firstParticleName+secondParticleName;
      _gFactorMap[i].push_back(currentGValue);
    }
  }

  std::map<int, std::vector<double> >::iterator itgFac;
  for (itgFac=_gFactorMap.begin(); itgFac!=_gFactorMap.end(); ++itgFac){
    std::vector<double> currentgVector=itgFac->second;
    std::shared_ptr<KPole> currentPole;
    if (_kMatrixParser->useBarrierFactors()) currentPole=std::shared_ptr<KPole>(new KPoleBarrier(currentgVector, _poleMasses.at(itgFac->first), _phpVecs, _kMatrixParser->orbitalMom()));
    else currentPole=std::shared_ptr<KPole>(new KPole(currentgVector, _poleMasses.at(itgFac->first)));
    _kPoles.push_back(currentPole);
  }

  int orderBg=_kMatrixParser->orderBg();
  if(orderBg<0) _kMatr=std::shared_ptr<KMatrixRel>(new KMatrixRel(_kPoles,_phpVecs ));
  else{
    bool withAdler=_kMatrixParser->useAdler();
    _kMatr=std::shared_ptr<KMatrixRel>(new KMatrixRelBg(_kPoles,_phpVecs, orderBg, withAdler));
    _kMatr->updateBgTerms(0, 0, 0,2.420422828649453);
    _kMatr->updateBgTerms(0, 0, 1,0.9134263692474104);
    _kMatr->updateBgTerms(0, 1, 1,0.9134546802616772);
    _kMatr->updateBgTerms(0, 0, 2,0.);
    _kMatr->updateBgTerms(0, 1, 2,0.);
    _kMatr->updateBgTerms(0, 2, 2,0.);
    // _kMatr->updateBgTerms(0, 0, 0, 0.0);
    // _kMatr->updateBgTerms(0, 0, 1, 0.0);
    // _kMatr->updateBgTerms(0, 1, 1, 0.0);
    if(orderBg>0){
      _kMatr->updateBgTerms(1, 0, 0, -0.15099);
      _kMatr->updateBgTerms(1, 0, 1, -0.038266);
      _kMatr->updateBgTerms(1, 1, 1, -0.0219);
    }
    if(orderBg>1){
      _kMatr->updateBgTerms(2, 0, 0, 0.00811);
      _kMatr->updateBgTerms(2, 0, 1, 0.0022596);
      _kMatr->updateBgTerms(2, 1, 1, 0.00085655);
    }

    if(withAdler){
      _kMatr->updates0Adler(_kMatrixParser->s0Adler());
      _kMatr->updatesnormAdler(_kMatrixParser->snormAdler());
    }
  }
  _tMatr=std::shared_ptr<TMatrixRel>(new TMatrixRel(_kMatr));


}

