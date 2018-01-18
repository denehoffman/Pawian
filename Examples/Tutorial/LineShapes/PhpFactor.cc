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
#include "Utils/PawianConstants.hh"
#include "PwaDynamics/PhaseSpaceIsobarLUT.hh"

#include "ErrLogger/ErrLogger.hh"

PhpFactor::PhpFactor(double mass1, double mass2, double massMax) :
  _theTFile(0)
  ,_mass1(mass1)
  ,_mass2(mass2)
  ,_massMax(massMax)
{
  _theTFile=new TFile("PhpFactor.root","recreate");

  PhaseSpaceIsobarLUT* phpBBUnstableRhoPi = new PhaseSpaceIsobarLUT(1.02, 0.135, "LUT/data/jollyj/mkuhlmann/PWA/LUT/LUTLoad/PhiPi2DimSqrtS.bin"); 

  _CMDefaultComplRealHist2= new TH2F("_CMDefaultComplRealHist2", "CM real default",601, -_massMax, _massMax, 601, -_massMax, _massMax);
  _CMDefaultComplRealHist2->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _CMDefaultComplRealHist2->GetXaxis()->SetTitle("Re(#sqrt{s})");
  _CMDefaultComplImagHist2= new TH2F("_CMDefaultComplImagHist2", "CM imag default",601, -_massMax, _massMax, 601, -_massMax, _massMax);
  _CMDefaultComplImagHist2->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _CMDefaultComplImagHist2->GetXaxis()->SetTitle("Re(#sqrt{s})");

  _CMPenningtonRealHist2= new TH2F("_CMPenningtonRealHist2", "CM real Pennington",601, -_massMax, _massMax, 601, -_massMax, _massMax);
  _CMPenningtonRealHist2->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _CMPenningtonRealHist2->GetXaxis()->SetTitle("Re(#sqrt{s})");
  _CMPenningtonImagHist2= new TH2F("_CMPenningtonImagHist2", "CM imag Pennington",601, -_massMax, _massMax, 601, -_massMax, _massMax);
  _CMPenningtonImagHist2->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _CMPenningtonImagHist2->GetXaxis()->SetTitle("Re(#sqrt{s})");

  _CMReidRealHist2= new TH2F("_CMReidRealHist2", "CM real Reid",601, -_massMax, _massMax, 601, -_massMax, _massMax);
  _CMReidRealHist2->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _CMReidRealHist2->GetXaxis()->SetTitle("Re(#sqrt{s})");
  _CMReidImagHist2= new TH2F("_CMReidImagHist2", "CM imag Reid",601, -_massMax, _massMax, 601, -_massMax, _massMax);
  _CMReidImagHist2->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _CMReidImagHist2->GetXaxis()->SetTitle("Re(#sqrt{s})");

  double massMaxSqr= _massMax* _massMax;  
  _CMAsnerRealHist= new TH1F("_CMAsnerRealHist", "CM real Asner",301, -massMaxSqr, massMaxSqr);
  _CMAsnerRealHist->GetYaxis()->SetTitle("Re(CM)");
  _CMAsnerRealHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");
  _CMAsnerImagHist= new TH1F("_CMAsnerImagHist", "CM imag Asner",301, -massMaxSqr, massMaxSqr);
  _CMAsnerImagHist->GetYaxis()->SetTitle("Im(CM)");
  _CMAsnerImagHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");

  _CMPenningtonRealHist= new TH1F("_CMPenningtonRealHist", "CM real Pennington",301, -massMaxSqr, massMaxSqr);
  _CMPenningtonRealHist->GetYaxis()->SetTitle("Re(CM)");
  _CMPenningtonRealHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");
  _CMPenningtonImagHist= new TH1F("_CMPenningtonImagHist", "CM imag Pennington",301, -massMaxSqr, massMaxSqr);
  _CMPenningtonImagHist->GetYaxis()->SetTitle("Im(CM)");
  _CMPenningtonImagHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");

  _CMReidRealHist= new TH1F("_CMReidRealHist", "CM real Reid",301, -massMaxSqr, massMaxSqr);
  _CMReidRealHist->GetYaxis()->SetTitle("Re(CM)");
  _CMReidRealHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");
  _CMReidImagHist= new TH1F("_CMReidImagHist", "CM imag Reid",301, -massMaxSqr, massMaxSqr);
  _CMReidImagHist->GetYaxis()->SetTitle("Im(CM)");
  _CMReidImagHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");

  _CMDefaultRealHist= new TH1F("_CMDefaultRealHist", "CM real default",301, -massMaxSqr, massMaxSqr);
  _CMDefaultRealHist->GetYaxis()->SetTitle("Re(CM)");
  _CMDefaultRealHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");
  _CMDefaultImagHist= new TH1F("_CMDefaultImagHist", "CM imag default",301, -massMaxSqr, massMaxSqr);
  _CMDefaultImagHist->GetYaxis()->SetTitle("Im(CM)");
  _CMDefaultImagHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");

  _CMDefaultComplRealHist= new TH1F("_CMDefaultComplRealHist", "CM real default",301, -massMaxSqr, massMaxSqr);
  _CMDefaultComplRealHist->GetYaxis()->SetTitle("Re(CM)");
  _CMDefaultComplRealHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");
  _CMDefaultComplImagHist= new TH1F("_CMDefaultComplImagHist", "#rho imag default",301, -massMaxSqr, massMaxSqr);
  _CMDefaultComplImagHist->GetYaxis()->SetTitle("Im(CM)");
  _CMDefaultComplImagHist->GetXaxis()->SetTitle("s[GeV^{2}/c^{4}]");




  _CMMAsnerRealHist= new TH1F("_CMMAsnerRealHist", "CM real Asner",301, -massMax, massMax);
  _CMMAsnerRealHist->GetYaxis()->SetTitle("Re(CM)");
  _CMMAsnerRealHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");
  _CMMAsnerImagHist= new TH1F("_CMMAsnerImagHist", "CM imag Asner",301, -massMax, massMax);
  _CMMAsnerImagHist->GetYaxis()->SetTitle("Im(CM)");
  _CMMAsnerImagHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");

  _CMMPenningtonRealHist= new TH1F("_CMMPenningtonRealHist", "CM real Pennington",301, -massMax, massMax);
  _CMMPenningtonRealHist->GetYaxis()->SetTitle("Re(CM)");
  _CMMPenningtonRealHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");
  _CMMPenningtonImagHist= new TH1F("_CMMPenningtonImagHist", "CM imag Pennington",301, -massMax, massMax);
  _CMMPenningtonImagHist->GetYaxis()->SetTitle("Im(CM)");
  _CMMPenningtonImagHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");

  _CMMReidRealHist= new TH1F("_CMMReidRealHist", "CM real Reid",301, -massMax, massMax);
  _CMMReidRealHist->GetYaxis()->SetTitle("Re(CM)");
  _CMMReidRealHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");
  _CMMReidImagHist= new TH1F("_CMMReidImagHist", "CM imag Reid",301, -massMax, massMax);
  _CMMReidImagHist->GetYaxis()->SetTitle("Im(CM)");
  _CMMReidImagHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");

  _CMMDefaultComplRealHist= new TH1F("_CMMDefaultComplRealHist", "CM real default",301, -massMax, massMax);
  _CMMDefaultComplRealHist->GetYaxis()->SetTitle("Re(CM)");
  _CMMDefaultComplRealHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");
  _CMMDefaultComplImagHist= new TH1F("_CMMDefaultComplImagHist", "CM imag default",301, -massMax, massMax);
  _CMMDefaultComplImagHist->GetYaxis()->SetTitle("Im(CM)");
  _CMMDefaultComplImagHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");

  _CMMBBUnstableRhoPiRealHist = new TH1F("_CMMBBUnstablePiPiRealHist", "CM real BB unstable phi pi",301, -massMax, massMax);
  _CMMBBUnstableRhoPiRealHist->GetYaxis()->SetTitle("Re(CM)");
  _CMMBBUnstableRhoPiRealHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");
  _CMMBBUnstableRhoPiImagHist = new TH1F("_CMMBUnstablePiPiImagHist", "CM imag BB unstable phi pi",301, -massMax, massMax);
  _CMMBBUnstableRhoPiImagHist->GetYaxis()->SetTitle("Im(CM)");
  _CMMBBUnstableRhoPiImagHist->GetXaxis()->SetTitle("m[GeV/c^{2}]");

  double stepSize=(2.*massMax)/301.;  

  for (double massIt=-_massMax+stepSize/2.; massIt<_massMax; massIt+=stepSize){
    complex<double> currentMass(massIt,0.);

    std::complex<double> currentCMFac = PawianConstants::i*PawianQFT::phaseSpaceFacDefault(currentMass, _mass1, _mass2);
    _CMMDefaultComplRealHist->Fill(massIt, currentCMFac.real());
    _CMMDefaultComplImagHist->Fill(massIt, currentCMFac.imag());

    complex<double> currentCMFacAsner = PawianQFT::ChewMandelstamAsner(massIt*massIt, _mass1, _mass2);
    _CMMAsnerRealHist->Fill(massIt, currentCMFacAsner.real());
    _CMMAsnerImagHist->Fill(massIt, currentCMFacAsner.imag());

    complex<double> currentMassSqr=currentMass*currentMass;
    complex<double> currentCMPenningtonFac = PawianQFT::ChewMandelstamPennington(currentMassSqr, _mass1, _mass2);
    _CMMPenningtonRealHist->Fill(massIt, currentCMPenningtonFac.real());
    _CMMPenningtonImagHist->Fill(massIt, currentCMPenningtonFac.imag());

    currentMass+=complex<double>(0.,1.e-10);
    complex<double> currentCMReidFac = PawianQFT::ChewMandelstamReid(currentMass*currentMass, _mass1, _mass2);
    _CMMReidRealHist->Fill(massIt, currentCMReidFac.real());
    _CMMReidImagHist->Fill(massIt, currentCMReidFac.imag());

    InfoMsg << "currentMass: " << massIt << endmsg;
    complex<double> currentCMBBUnstableRhoPiFac = phpBBUnstableRhoPi->ChewM(massIt);
    _CMMBBUnstableRhoPiRealHist->Fill(massIt, currentCMBBUnstableRhoPiFac.real());
    _CMMBBUnstableRhoPiImagHist->Fill(massIt, currentCMBBUnstableRhoPiFac.imag());    
  }

  stepSize=2.*massMax/601.;  
  for (double massRealIt=-_massMax+stepSize/2.; massRealIt<_massMax; massRealIt+=stepSize){
    for (double massImagIt=-_massMax+stepSize/2.; massImagIt<_massMax; massImagIt+=stepSize){
      complex<double> currentMass(massRealIt, massImagIt);

      std::complex<double> currentCMFac = PawianConstants::i*PawianQFT::phaseSpaceFacDefault(currentMass, _mass1, _mass2);
      _CMDefaultComplRealHist2->Fill(massRealIt, massImagIt, currentCMFac.real());
      _CMDefaultComplImagHist2->Fill(massRealIt, massImagIt, currentCMFac.imag());

      complex<double> currentMassSqr=currentMass*currentMass;
      complex<double> currentCMPenningtonFac = PawianQFT::ChewMandelstamPennington(currentMassSqr, _mass1, _mass2);
      _CMPenningtonRealHist2->Fill(massRealIt, massImagIt, currentCMPenningtonFac.real());
      _CMPenningtonImagHist2->Fill(massRealIt, massImagIt, currentCMPenningtonFac.imag());

      complex<double> currentCMReidFac = PawianQFT::ChewMandelstamReid(currentMassSqr, _mass1, _mass2);
      _CMReidRealHist2->Fill(massRealIt, massImagIt, currentCMReidFac.real());
      _CMReidImagHist2->Fill(massRealIt, massImagIt, currentCMReidFac.imag());
    }
  }

  stepSize=2.*massMaxSqr/301.;
  for (double massSqrIt=-massMaxSqr+stepSize/2.; massSqrIt < massMaxSqr+stepSize/2.; massSqrIt+=stepSize){
    complex<double> currentCMFac = PawianQFT::ChewMandelstamAsner(massSqrIt, _mass1, _mass2);
    if(currentCMFac.real()==currentCMFac.real() && currentCMFac.imag()==currentCMFac.imag()){ 
    _CMAsnerRealHist->Fill(massSqrIt, currentCMFac.real());
    _CMAsnerImagHist->Fill(massSqrIt, currentCMFac.imag());
    }

    complex<double> currentCMDefaultFac = PawianConstants::i*PawianQFT::phaseSpaceFacDefault(massSqrIt, _mass1, _mass2);
    if(currentCMDefaultFac.real()==currentCMDefaultFac.real() && currentCMDefaultFac.imag()==currentCMDefaultFac.imag() && fabs(massSqrIt)>1.e-3){
      _CMDefaultRealHist->Fill(massSqrIt, currentCMDefaultFac.real());
      _CMDefaultImagHist->Fill(massSqrIt, currentCMDefaultFac.imag());
    }

    complex<double> massSqrItCompl(massSqrIt, 0.);
    complex<double> currentCMDefaultComplFac =  PawianConstants::i*PawianQFT::phaseSpaceFacDefault(massSqrItCompl, _mass1, _mass2);
  
    _CMDefaultComplRealHist->Fill(massSqrItCompl.real(), currentCMDefaultComplFac.real());
    _CMDefaultComplImagHist->Fill(massSqrItCompl.real(), currentCMDefaultComplFac.imag());


    complex<double> currentCMPenningtonFac = PawianQFT::ChewMandelstamPennington(massSqrItCompl, _mass1, _mass2);
    _CMPenningtonRealHist->Fill(massSqrItCompl.real(), currentCMPenningtonFac.real());
    _CMPenningtonImagHist->Fill(massSqrItCompl.real(), currentCMPenningtonFac.imag());

    complex<double> massSqrItComplCM=massSqrItCompl+complex<double>(0., 0.0000001); 
    complex<double> currentReidFac = PawianQFT::ChewMandelstamReid(massSqrItComplCM, _mass1, _mass2);
    _CMReidRealHist->Fill(massSqrItCompl.real(), currentReidFac.real());
    _CMReidImagHist->Fill(massSqrItCompl.real(), currentReidFac.imag());

  }
  delete phpBBUnstableRhoPi; 
}

PhpFactor::~PhpFactor()
{
   _theTFile->Write();
   _theTFile->Close();
}

