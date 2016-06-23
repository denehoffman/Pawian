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
#include <boost/multi_array.hpp>
#include "Examples/Tutorial/LineShapes/TMatrixPiPiSWaveSimple4piPhp.hh"
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
#include "PwaDynamics/KMatrixPiPiSsimple4piPhp.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include "PwaDynamics/PhaseSpace4Pi.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/PdtParser.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

MatrixPiPiSWaveSimple4piPhp::MatrixPiPiSWaveSimple4piPhp(int numStepsForSheetScan, std::vector<double> energyPlaneBorders) :
  _noOfSteps(1000)
  ,_stepSize(0.)
  ,_massMin(.3)
  ,_massMax(1.5)
  ,_noOfChannels(5)
  ,_theTFile(0)
{ 
  init();
  std::string rootFileName="./MatrixPiPiSWaveSimple4piPhp.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");

  _massMax+=0.4;
  _stepSize=(_massMax-_massMin)/_noOfSteps;

  std::string ampRealKey="AmpReal";
  std::string ampImagKey="AmpImag";
  std::string argandKey="Argand";
  std::string phaseKey="Phase";
  std::string inelasticityKey="Elasticity";

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
	   << "\nenergyPlaneBorders 0:" << energyPlaneBorders[0] 
	   << "\nenergyPlaneBorders 1:" << energyPlaneBorders[1]
	   << "\nenergyPlaneBorders 2:" << energyPlaneBorders[2] 
	   << "\nenergyPlaneBorders 3:" << energyPlaneBorders[3] 
  	   << endmsg;
  std::vector<std::shared_ptr<AbsPhaseSpace> > phpVecs=_tMatr->kMatrix()->phaseSpaceVec();

  for (double mass=_massMin+_stepSize/0.5; mass<_massMax; mass+=_stepSize){
    Vector4<double> mass4Vec(mass, 0.,0.,0.);
    _tMatr->evalMatrix(mass);
    for(unsigned int i=0; i<_gFactorNames.size(); ++i){
      complex<double> currentRho=phpVecs.at(i)->factor(mass);

      _AmpRealH1Vec.at(i)->Fill(mass, sqrt(currentRho.real())*(*_tMatr)(i,i).real());
      _AmpImagH1Vec.at(i)->Fill(mass, sqrt(currentRho.real())*(*_tMatr)(i,i).imag());
      
      _ArgandH2Vec.at(i)->Fill( currentRho.real()*(*_tMatr)(i,i).real(), currentRho.real()*(*_tMatr)(i,i).imag());
      double currentphase=360.*atan2((*_tMatr)(i,i).imag(),(*_tMatr)(i,i).real()) / 3.1415;
      _PhaseH2Vec.at(i)->Fill(mass, currentphase);
      double sqrtFactor=(*_tMatr)(i,i).real()*(*_tMatr)(i,i).real()+((*_tMatr)(i,i).imag()-0.5)*((*_tMatr)(i,i).imag()-0.5);
      double currentElasticity=2.*sqrt(sqrtFactor);
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

  RiemannSheetAnalyzer(_noOfChannels, _tMatr,
		       std::complex<double>(energyPlaneBorders[0], energyPlaneBorders[1]),
		       std::complex<double>(energyPlaneBorders[2], energyPlaneBorders[3]),
		       numStepsForSheetScan);
}

MatrixPiPiSWaveSimple4piPhp::~MatrixPiPiSWaveSimple4piPhp()
{
  _theTFile->Write();
  _theTFile->Close();
}


void MatrixPiPiSWaveSimple4piPhp::init(){
   PdtParser pdtParser;
   std::string theSourcePath=getenv("TOP_DIR");
   std::string pdtFileRelPath="/Particle/pdtNew.table";
   std::string pdtFile(theSourcePath+pdtFileRelPath);
   _particleTable = new ParticleTable;

   if (!pdtParser.parse(pdtFile, *_particleTable)) {
      Alert << "can not parse particle table " << pdtFile << endmsg;
      exit(1);
   }

   _poleNames.push_back("pole1");
   _poleNames.push_back("pole2");
   _poleNames.push_back("pole3");
   _poleNames.push_back("pole4");
   _poleNames.push_back("pole5");


    _gFactorNames.push_back("g1"); 
    _gFactorNames.push_back("g2");     
    _gFactorNames.push_back("g3");
    _gFactorNames.push_back("g4");     
    _gFactorNames.push_back("g5");

    _kMatr=std::shared_ptr<KMatrixBase>(new KMatrixPiPiSsimple4piPhp()); 
    _tMatr=std::shared_ptr<TMatrixRel>(new TMatrixRel(_kMatr));
}

