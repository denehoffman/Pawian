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
#include <complex>
#include "Examples/Tutorial/LineShapes/PhpFactor.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include "qft++Extension/PawianUtils.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

PhpFactor::PhpFactor(double mass1, double mass2, double massMax) :
  _theTFile(0)
  ,_mass1(mass1)
  ,_mass2(mass2)
  ,_massMax(massMax)
{
  complex<double> i(0.,1.);
  _theTFile=new TFile("PhpFactor.root","recreate");
 
  _phpRealHist= new TH1F("_phpRealHist", "#rho real",301, -_massMax, _massMax);
  _phpImagHist= new TH1F("_phpImagHist", "#rho imag",301, -_massMax, _massMax);

  _phpDefaultComplRealHist2= new TH2F("_phpDefaultComplRealHist2", "#rho real default",601, -_massMax, _massMax, 601, -_massMax, _massMax);
  _phpDefaultComplRealHist2->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _phpDefaultComplRealHist2->GetXaxis()->SetTitle("Re(#sqrt{s})");
  _phpDefaultComplImagHist2= new TH2F("_phpDefaultComplImagHist2", "#rho imag default",601, -_massMax, _massMax, 601, -_massMax, _massMax);
  _phpDefaultComplImagHist2->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _phpDefaultComplImagHist2->GetXaxis()->SetTitle("Re(#sqrt{s})");

  _phpPenningtonRealHist2= new TH2F("_phpPenningtonRealHist2", "#rho real Pennington",601, -_massMax, _massMax, 601, -_massMax, _massMax);
  _phpPenningtonRealHist2->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _phpPenningtonRealHist2->GetXaxis()->SetTitle("Re(#sqrt{s})");
  _phpPenningtonImagHist2= new TH2F("_phpPenningtonImagHist2", "#rho imag Pennington",601, -_massMax, _massMax, 601, -_massMax, _massMax);
  _phpPenningtonImagHist2->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _phpPenningtonImagHist2->GetXaxis()->SetTitle("Re(#sqrt{s})");

  _phpReidRealHist2= new TH2F("_phpReidRealHist2", "#rho real Reid",601, -_massMax, _massMax, 601, -_massMax, _massMax);
  _phpReidRealHist2->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _phpReidRealHist2->GetXaxis()->SetTitle("Re(#sqrt{s})");
  _phpReidImagHist2= new TH2F("_phpReidImagHist2", "#rho imag Reid",601, -_massMax, _massMax, 601, -_massMax, _massMax);
  _phpReidImagHist2->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _phpReidImagHist2->GetXaxis()->SetTitle("Re(#sqrt{s})");

  double massMaxSqr= _massMax* _massMax;  
  _phpAsnerRealHist= new TH1F("_phpAsnerRealHist", "#rho real Asner",301, -massMaxSqr, massMaxSqr);
  _phpAsnerRealHist->GetYaxis()->SetTitle("i #rho");
  _phpAsnerRealHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");
  _phpAsnerImagHist= new TH1F("_phpAsnerImagHist", "#rho imag Asner",301, -massMaxSqr, massMaxSqr);
  _phpAsnerImagHist->GetYaxis()->SetTitle("i #rho");
  _phpAsnerImagHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");

  _phpPenningtonRealHist= new TH1F("_phpPenningtonRealHist", "#rho real Pennington",301, -massMaxSqr, massMaxSqr);
  _phpPenningtonRealHist->GetYaxis()->SetTitle("i #rho");
  _phpPenningtonRealHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");
  _phpPenningtonImagHist= new TH1F("_phpPenningtonImagHist", "#rho imag Pennington",301, -massMaxSqr, massMaxSqr);
  _phpPenningtonImagHist->GetYaxis()->SetTitle("i #rho");
  _phpPenningtonImagHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");

  _phpReidRealHist= new TH1F("_phpReidRealHist", "#rho real Reid",301, -massMaxSqr, massMaxSqr);
  _phpReidRealHist->GetYaxis()->SetTitle("i #rho");
  _phpReidRealHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");
  _phpReidImagHist= new TH1F("_phpReidImagHist", "#rho imag Reid",301, -massMaxSqr, massMaxSqr);
  _phpReidImagHist->GetYaxis()->SetTitle("i #rho");
  _phpReidImagHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");

  _phpDefaultRealHist= new TH1F("_phpDefaultRealHist", "#rho real default",301, -massMaxSqr, massMaxSqr);
  _phpDefaultImagHist= new TH1F("_phpDefaultImagHist", "#rho imag default",301, -massMaxSqr, massMaxSqr);
  _phpDefaultComplRealHist= new TH1F("_phpDefaultComplRealHist", "#rho real default",301, -massMaxSqr, massMaxSqr);
  _phpDefaultComplRealHist->GetYaxis()->SetTitle("i #rho");
  _phpDefaultComplRealHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");
  _phpDefaultComplImagHist= new TH1F("_phpDefaultComplImagHist", "#rho imag default",301, -massMaxSqr, massMaxSqr);
  _phpDefaultComplImagHist->GetYaxis()->SetTitle("i #rho");
  _phpDefaultComplImagHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");




  _phpMAsnerRealHist= new TH1F("_phpMAsnerRealHist", "#rho real Asner",301, -massMax, massMax);
  _phpMAsnerRealHist->GetYaxis()->SetTitle("i #rho");
  _phpMAsnerRealHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");
  _phpMAsnerImagHist= new TH1F("_phpMAsnerImagHist", "#rho imag Asner",301, -massMax, massMax);
  _phpMAsnerImagHist->GetYaxis()->SetTitle("i #rho");
  _phpMAsnerImagHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");

  _phpMPenningtonRealHist= new TH1F("_phpMPenningtonRealHist", "#rho real Pennington",301, -massMax, massMax);
  _phpMPenningtonRealHist->GetYaxis()->SetTitle("i #rho");
  _phpMPenningtonRealHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");
  _phpMPenningtonImagHist= new TH1F("_phpMPenningtonImagHist", "#rho imag Pennington",301, -massMax, massMax);
  _phpMPenningtonImagHist->GetYaxis()->SetTitle("i #rho");
  _phpMPenningtonImagHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");

  _phpMReidRealHist= new TH1F("_phpMReidRealHist", "#rho real Reid",301, -massMax, massMax);
  _phpMReidRealHist->GetYaxis()->SetTitle("i #rho");
  _phpMReidRealHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");
  _phpMReidImagHist= new TH1F("_phpMReidImagHist", "#rho imag Reid",301, -massMax, massMax);
  _phpMReidImagHist->GetYaxis()->SetTitle("i #rho");
  _phpMReidImagHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");

  _phpMDefaultComplRealHist= new TH1F("_phpMDefaultComplRealHist", "#rho real default",301, -massMax, massMax);
  _phpMDefaultComplRealHist->GetYaxis()->SetTitle("i #rho");
  _phpMDefaultComplRealHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");
  _phpMDefaultComplImagHist= new TH1F("_phpMDefaultComplImagHist", "#rho imag default",301, -massMax, massMax);
  _phpMDefaultComplImagHist->GetYaxis()->SetTitle("i #rho");
  _phpMDefaultComplImagHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");

  double stepSize=(2.*massMax)/301.;  

  for (double massIt=-_massMax+stepSize/2.; massIt<_massMax; massIt+=stepSize){
    complex<double> currentMass(massIt,0.);
    //    if (norm(currentMass) < 1.e-5) continue;
  
    //    std::complex<double> currentFac = phaseSpaceFac(currentMass, _mass1, _mass2);
    std::complex<double> currentFac = PawianQFT::phaseSpaceFacDefault(currentMass, _mass1, _mass2);
    std::complex<double> currentFaci = i*currentFac;
    double currentImag = currentFaci.imag();
    if (currentImag>1000.) currentImag=1000.;
    _phpRealHist->Fill(massIt, currentFaci.real());
    _phpImagHist->Fill(massIt, currentImag);

    _phpMDefaultComplRealHist->Fill(massIt,currentFaci.real());
    _phpMDefaultComplImagHist->Fill(massIt,currentImag);

    complex<double> currentFacAsner = PawianQFT::phaseSpaceFacAsner(massIt*massIt, _mass1, _mass2);
    _phpMAsnerRealHist->Fill(massIt, currentFacAsner.real());
    _phpMAsnerImagHist->Fill(massIt, currentFacAsner.imag());

    complex<double> currentMassSqr=currentMass*currentMass;
    complex<double> currentPenningtonFac = PawianQFT::phaseSpaceFacPennington(currentMassSqr, _mass1, _mass2);
    _phpMPenningtonRealHist->Fill(massIt, currentPenningtonFac.real());
    _phpMPenningtonImagHist->Fill(massIt, currentPenningtonFac.imag());

    complex<double> currentReidFac = i*PawianQFT::phaseSpaceFacReid(massIt, _mass1, _mass2);
    //    if ( fabs(currentReidFac.real()) <1000. &&  fabs(currentReidFac.imag()) <1000.){ 
      _phpMReidRealHist->Fill(massIt, currentReidFac.real());
      _phpMReidImagHist->Fill(massIt, currentReidFac.imag());    
      //    }
  }

  stepSize=2.*massMax/601.;  
  for (double massRealIt=-_massMax+stepSize/2.; massRealIt<_massMax; massRealIt+=stepSize){
    for (double massImagIt=-_massMax+stepSize/2.; massImagIt<_massMax; massImagIt+=stepSize){
      complex<double> currentMass(massRealIt, massImagIt);
      //      std::complex<double> currentFac = i*phaseSpaceFac(currentMass, _mass1, _mass2);
      std::complex<double> currentFac = i*PawianQFT::phaseSpaceFacDefault(currentMass, _mass1, _mass2);
      _phpDefaultComplRealHist2->Fill(massRealIt, massImagIt, currentFac.real());
      _phpDefaultComplImagHist2->Fill(massRealIt, massImagIt, currentFac.imag());

      complex<double> currentMassSqr=currentMass*currentMass;
      complex<double> currentPenningtonFac = PawianQFT::phaseSpaceFacPennington(currentMassSqr, _mass1, _mass2);
      _phpPenningtonRealHist2->Fill(massRealIt, massImagIt, currentPenningtonFac.real());
      _phpPenningtonImagHist2->Fill(massRealIt, massImagIt, currentPenningtonFac.imag());

      complex<double> currentReidFac = i*PawianQFT::phaseSpaceFacReid(currentMass, _mass1, _mass2);
      _phpReidRealHist2->Fill(massRealIt, massImagIt, currentReidFac.real());
      _phpReidImagHist2->Fill(massRealIt, massImagIt, currentReidFac.imag());
    }
  }

  stepSize=2.*massMaxSqr/301.;
  for (double massSqrIt=-massMaxSqr+stepSize/2.; massSqrIt < massMaxSqr+stepSize/2.; massSqrIt+=stepSize){
    complex<double> currentFac = PawianQFT::phaseSpaceFacAsner(massSqrIt, _mass1, _mass2);
    // std::cout << "massSqrIt: " << massSqrIt << std::endl; 
    // std::cout << "currentFac.real(): " << currentFac.real() << std::endl;
    // std::cout << "currentFac.imag(): " << currentFac.imag() << std::endl;
    if(currentFac.real()==currentFac.real() && currentFac.imag()==currentFac.imag()){ 
    _phpAsnerRealHist->Fill(massSqrIt, currentFac.real());
    _phpAsnerImagHist->Fill(massSqrIt, currentFac.imag());
    }

    //    complex<double> currentDefaultFac = i*phaseSpaceFac(massSqrIt, _mass1, _mass2);
    complex<double> currentDefaultFac = i*PawianQFT::phaseSpaceFacDefault(massSqrIt, _mass1, _mass2);
    if(currentDefaultFac.real()==currentDefaultFac.real() && currentDefaultFac.imag()==currentDefaultFac.imag() && fabs(massSqrIt)>1.e-3){
      _phpDefaultRealHist->Fill(massSqrIt, currentDefaultFac.real());
      _phpDefaultImagHist->Fill(massSqrIt, currentDefaultFac.imag());
    }

    complex<double> massSqrItCompl(massSqrIt, 0.);
    //    complex<double> currentDefaultComplFac = i*phaseSpaceFac(massSqrItCompl, _mass1, _mass2);
    complex<double> currentDefaultComplFac = i*PawianQFT::phaseSpaceFacDefault(massSqrItCompl, _mass1, _mass2);
    if(currentDefaultComplFac.real()==currentDefaultComplFac.real() && currentDefaultComplFac.imag()==currentDefaultComplFac.imag()){
      double currentDefaultComplFacReal=currentDefaultComplFac.real();
      if (currentDefaultComplFacReal<-1000.) currentDefaultComplFacReal=-1000.;
      else if (currentDefaultComplFacReal>1000.) currentDefaultComplFacReal=1000.;

      double currentDefaultComplFacImag=currentDefaultComplFac.imag();
      if (currentDefaultComplFacImag<-1000.) currentDefaultComplFacImag=-1000.;
      else if (currentDefaultComplFacImag>1000.) currentDefaultComplFacImag=1000.;
  
      _phpDefaultComplRealHist->Fill(massSqrItCompl.real(), currentDefaultComplFacReal);
      _phpDefaultComplImagHist->Fill(massSqrItCompl.real(), currentDefaultComplFacImag);
    }

    complex<double> currentPenningtonFac = PawianQFT::phaseSpaceFacPennington(massSqrItCompl, _mass1, _mass2);
    if(currentPenningtonFac.real()==currentPenningtonFac.real() && currentPenningtonFac.imag()==currentPenningtonFac.imag() && fabs(currentPenningtonFac.real())>1.e-3){
      _phpPenningtonRealHist->Fill(massSqrItCompl.real(), currentPenningtonFac.real());
      _phpPenningtonImagHist->Fill(massSqrItCompl.real(), currentPenningtonFac.imag());
    }

    complex<double> massSqrItComplCM=massSqrItCompl+complex<double>(0., 0.0000001); 
    complex<double> currentReidFac = PawianQFT::ChewMandelstamReid(massSqrItComplCM, _mass1, _mass2);
    _phpReidRealHist->Fill(massSqrItCompl.real(), currentReidFac.real());
    _phpReidImagHist->Fill(massSqrItCompl.real(), currentReidFac.imag());        
  }  
}

PhpFactor::~PhpFactor()
{
   _theTFile->Write();
   _theTFile->Close();
}

