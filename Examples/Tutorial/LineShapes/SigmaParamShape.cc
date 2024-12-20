//************************************************************************//
//									  //
//  Copyright 2024 Bertram Kopf (bertram@ep1.rub.de)			  //
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
#include "Examples/Tutorial/LineShapes/SigmaParamShape.hh"
#include "PwaDynamics/SigmaParameterization.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include "Utils/PawianConstants.hh"
#include "qft++Extension/PawianUtils.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

SigmaParamShape::SigmaParamShape() :
  _theTFile(0)
  ,_graphPhase(new TGraph())
  ,_graphIntensity(new TGraph())
  ,_graphArgand(new TGraph())
{
  std::string rootFileName="./SigmaShape.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");

  _graphPhase->SetName("grPhase");
  _graphPhase->SetTitle("#sigma phase #pi #pi");

  _graphIntensity->SetName("grIntensity");
  _graphIntensity->SetTitle("#sigma intensity #pi #pi");

  _graphArgand->SetName("grArgand");
  _graphArgand->SetTitle("#sigma Argand #pi #pi");
  
  int size=100;
  double massMin=.3;
  double massMax=1.2;

  double stepSize=(massMax-massMin)/size;

  //  _histPhase= new TH1F("_histPhase","hist phase",size+1, massMin, massMax);


   SigmaParameterization theSigmaParam; 
   int counter=0;  
   for (double mass=massMin; mass<massMax; mass+=stepSize){
     ++counter;
     complex<double> currentAmp=theSigmaParam.calc(mass);
     //     _histPhase->Fill(std::arg(currentAmp)*180./TMath::Pi());
     _graphPhase->SetPoint(counter, mass, std::arg(currentAmp)*180./TMath::Pi());

     //     std::complex<double> rhopipi=PawianQFT::phaseSpaceFacDefault(mass,PawianConstants::mPi, PawianConstants::mPi);
     //     std::complex<double> breapUpMom=PawianQFT::breakupMomQDefault(mass, PawianConstants::mPi, PawianConstants::mPi);
     _graphIntensity->SetPoint(counter, mass, std::norm(currentAmp));

     std::complex<double> Tmat=theSigmaParam.calcT(mass); 
     //_graphArgand->SetPoint(counter, sqrt(rhopipi.real())*Tmat.real(), sqrt(rhopipi.real())*Tmat.imag());
     _graphArgand->SetPoint(counter, Tmat.real(), Tmat.imag());

   }
   
}

SigmaParamShape::~SigmaParamShape()
{
  _graphPhase->Write();
  _graphIntensity->Write();
  _graphArgand->Write();
  _theTFile->Write();
  _theTFile->Close();
}

