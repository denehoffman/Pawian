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
#include "TGraph.h"
#include "TMath.h"
#include "Utils/PawianConstants.hh"
#include "PwaDynamics/PhaseSpaceIsobarLUT.hh"

#include "ErrLogger/ErrLogger.hh"

LUTPlot::LUTPlot(std::string LUTFilePath, double realGranularity, double imagGranularity) :
  _theTFile(new TFile("LUTPlot.root","recreate"))
  ,_graphReal(new TGraph())
  ,_graphImag(new TGraph())
  ,_graph2DReal(new TGraph2D())
  ,_graph2DImag(new TGraph2D())
  ,_LUTFilePath(LUTFilePath)
  ,_realGranularity(realGranularity)
  ,_imagGranularity(imagGranularity)
{
  PhaseSpaceIsobarLUT* theFactor = new PhaseSpaceIsobarLUT(0.0, 0.0, _LUTFilePath.c_str()); 
  const std::complex<double> upperLimit = theFactor->sHigh();
  const std::complex<double> lowerLimit = theFactor->sLow();

  const int stepCountReal = abs((int)((upperLimit.real()-lowerLimit.real())/realGranularity));
  const int stepCountImag = abs((int)((lowerLimit.imag()-upperLimit.imag())/imagGranularity));

  _graphReal->SetTitle("BBUnstable real part");
  _graphImag->SetTitle("BBUnstable imaginary part");

  _graph2DReal->SetTitle("BBUnstable real part 2D");
  _graph2DImag->SetTitle("BBUnstable imaginary part 2D");
  
  for(int i = 1; i < stepCountReal; i++){
    complex<double> currentMass(lowerLimit.real()+i*realGranularity,0.0);
    // if(std::abs(currentMass.real())<2.*0.135) continue;
    complex<double> currentCMBBUnstableFac = theFactor->ChewM(currentMass);
    _graphReal->SetPoint(i, currentMass.real(), currentCMBBUnstableFac.real());
  }

  for(int i = 1; i < stepCountReal; i++){
    complex<double> currentMass(lowerLimit.real()+i*realGranularity,0.0);
    // if(std::abs(currentMass.real()) < 2.*0.135) continue;
    complex<double> currentCMBBUnstableFac = theFactor->ChewM(currentMass);
    _graphImag->SetPoint(i, currentMass.real(), currentCMBBUnstableFac.imag());
  }

  int counter=1;
  for(int i = 1; i < stepCountReal; i++){
    for(int j = 1; j < stepCountImag; j++){
      complex<double> currentMass(lowerLimit.real()+i*realGranularity,upperLimit.imag()+j*realGranularity);
      // if(std::abs(currentMass.real())<2.*0.135) continue;
      // if(std::abs(currentMass.imag())<0.005) continue;
      complex<double> currentCMBBUnstableFac = theFactor->ChewM(currentMass);
      _graph2DReal->SetPoint(counter,currentMass.real(), currentMass.imag(), currentCMBBUnstableFac.real());
      ++counter;
    }
  }

  counter=1;
  for(int i = 1; i < stepCountReal; i++){
    for(int j = 1; j < stepCountImag; j++){
      complex<double> currentMass(lowerLimit.real()+i*realGranularity,upperLimit.imag()+j*realGranularity);
      // if(std::abs(currentMass.real())<2.*0.135) continue;
      // if(std::abs(currentMass.imag())<0.005) continue;
      complex<double> currentCMBBUnstableFac = theFactor->ChewM(currentMass);
      _graph2DImag->SetPoint(counter, currentMass.real(), currentMass.imag(), currentCMBBUnstableFac.imag());
      ++counter;
    }
  }
  delete theFactor; 
}

LUTPlot::~LUTPlot()
{
  _graphReal->Write();
  _graphImag->Write();
  _graph2DReal->Write();
  _graph2DImag->Write();
   _theTFile->Write();
   _theTFile->Close();
}

