//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//		   Markus Kuhlmann (mkuhlmann@ep1.rub.de)		  //
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
#include "Examples/Tutorial/LineShapes/LUTPlot.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include "qft++Extension/PawianUtils.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "Utils/PawianConstants.hh"
#include "PwaDynamics/PhaseSpaceIsobarLUT.hh"

#include "ErrLogger/ErrLogger.hh"

LUTPlot::LUTPlot(std::string LUTFilePath, double realGranularity, double imagGranularity) :
  _theTFile(0)
  ,_LUTFilePath(LUTFilePath)
  ,_realGranularity(realGranularity)
  ,_imagGranularity(imagGranularity)
{
  _theTFile=new TFile("LUTPlot.root","recreate");

  PhaseSpaceIsobarLUT* theFactor = new PhaseSpaceIsobarLUT(0.0, 0.0, _LUTFilePath.c_str()); 
  const std::complex<double> upperLimit = theFactor->sHigh();
  const std::complex<double> lowerLimit = theFactor->sLow();

  const int stepCountReal = abs((int)((upperLimit.real()-lowerLimit.real())/realGranularity));
  const int stepCountImag = abs((int)((lowerLimit.imag()-upperLimit.imag())/imagGranularity));

  _RealPart2D = new TH2F("_RealPart2D", "BBUnstable 2D real part",stepCountReal-1, lowerLimit.real(), upperLimit.real()-0.01, stepCountImag-1, upperLimit.imag()+0.01, lowerLimit.imag()-0.01);
  _RealPart2D->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _RealPart2D->GetXaxis()->SetTitle("Re(#sqrt{s})");
  _RealPart2D->GetZaxis()->SetTitle("Re(CM)");

  _ImagPart2D = new TH2F("_ImagPart2D", "BBUnstable 2D imaginary part",stepCountReal-1, lowerLimit.real(), upperLimit.real()-0.01, stepCountImag-1, upperLimit.imag()+0.01, lowerLimit.imag()-0.01);
  _ImagPart2D->GetYaxis()->SetTitle("Im(#sqrt{s})");
  _ImagPart2D->GetXaxis()->SetTitle("Re(#sqrt{s})");
  _ImagPart2D->GetZaxis()->SetTitle("Im(CM)");

  _RealPart = new TH1F("_RealPart", "BBUnstable real part",stepCountReal-1, lowerLimit.real(), upperLimit.real()-0.01);
  _RealPart->GetYaxis()->SetTitle("Re(CM)");
  _RealPart->GetXaxis()->SetTitle("Re(#sqrt{s})");

  _ImagPart = new TH1F("_ImagPart", "BBUnstable imaginary part",stepCountReal-1, lowerLimit.real(), upperLimit.real()-0.01);
  _ImagPart->GetYaxis()->SetTitle("Im(CM)");
  _ImagPart->GetXaxis()->SetTitle("Re(#sqrt{s})");



  for(int i = 1; i < stepCountReal; i++){
    complex<double> currentMass(lowerLimit.real()+i*realGranularity,0.0);
    complex<double> currentCMBBUnstableFac = theFactor->ChewM(currentMass);
    _RealPart->Fill(currentMass.real(), currentCMBBUnstableFac.real());
  }

  for(int i = 1; i < stepCountReal; i++){
    complex<double> currentMass(lowerLimit.real()+i*realGranularity,0.0);
    complex<double> currentCMBBUnstableFac = theFactor->ChewM(currentMass);
    _ImagPart->Fill(currentMass.real(), currentCMBBUnstableFac.imag());
  }

  for(int i = 1; i < stepCountReal; i++){
    for(int j = 1; j < stepCountImag; j++){
      complex<double> currentMass(lowerLimit.real()+i*realGranularity,upperLimit.imag()+j*imagGranularity);
      complex<double> currentCMBBUnstableFac = theFactor->ChewM(currentMass);
      _RealPart2D->Fill(currentMass.real(), currentMass.imag(), currentCMBBUnstableFac.real());
    }
  }

  for(int i = 1; i < stepCountReal; i++){
    for(int j = 1; j < stepCountImag; j++){
      complex<double> currentMass(lowerLimit.real()+i*realGranularity,upperLimit.imag()+j*imagGranularity);
      complex<double> currentCMBBUnstableFac = theFactor->ChewM(currentMass);
      _ImagPart2D->Fill(currentMass.real(), currentMass.imag(), -currentCMBBUnstableFac.imag());
    }
  }
  delete theFactor; 
}

LUTPlot::~LUTPlot()
{
   _theTFile->Write();
   _theTFile->Close();
}

