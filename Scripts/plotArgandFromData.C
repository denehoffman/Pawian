//************************************************************************//
//                                                                        //
//  Copyright 2021 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 - Ruhr-Universität Bochum                              //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TKey.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TExec.h"

#include <algorithm> 
#include <vector>
#include <complex>

TGraph* argandGraph;
std::vector<float> massResultVec;

TCanvas* canvas = new TCanvas("Argand", "Argand", 20, 20, 900,900);


void drawtext() {
  Int_t i, n;
  Double_t x, y;
  TLatex *l;
  
  n = argandGraph->GetN();
  double oldX = 100.;
  double oldY = 100.;
  for (i=1; i<n; i++) {
    argandGraph->GetPoint(i, x, y);
    double currentDist=sqrt((oldX-x) * (oldX-x) + (oldY-y) * (oldY-y));
    if( currentDist > 0.15) {
      oldX = x;
      oldY = y;    
      l = new TLatex(x, y+0.05, Form("%4.2f", massResultVec.at(i)));
      l->SetTextSize(0.025);
      l->SetTextFont(42);
      l->SetTextAlign(21);
      l->Paint();
    }
  }
}

void plotArgandFromData(std::string iMagFileName="/data/jollyj/bertram/pwaData/ScatteringNew/I0/Swave/eta00_Madrid0.01sys_Ochs2013---SmallerErrors.out", std::string iPhaseFileName="/data/jollyj/bertram/pwaData/ScatteringNew/I0/Swave/delta00_Madrid0.01sys_Ochs2013---SmallerErrors.out", float matchAccuracy=0.01) {
  std::ifstream iMagStream(iMagFileName.c_str());
  if (!iMagStream) {
    cout << "can not open " << iMagFileName << endl;  // << endmsg;
    exit(1);
  }
  std::ifstream iPhaseStream(iPhaseFileName.c_str());
  if (!iPhaseStream) {
    cout << "can not open " << iPhaseFileName << endl;  // << endmsg;                                    
    exit(1);
  }
  std::vector<float> phaseVec;
  std::vector<float> massPhaseVec;
  std::vector<float> etaVec;
  std::vector<float> massEtaVec;

  while (!iMagStream.eof()) {
    double m, eta, etaErr;
    iMagStream >> m >> eta >> etaErr;
    cout << "m: " << m << " eta: " << eta << " etaErr: " << etaErr << endl;
    massEtaVec.push_back(m);
    etaVec.push_back(eta);
  }

  while (!iPhaseStream.eof()) {
    double m, phase, phaseErr;
    iPhaseStream >> m >> phase >> phaseErr;
    cout << "m: " << m << " phase: " << phase << " phaseErr: " << phaseErr << endl;
    massPhaseVec.push_back(m);
    phaseVec.push_back(phase);
  }
  
  std::vector<float> phaseResultVec;
  std::vector<float> etaResultVec;

  for(unsigned int i=0; i<massEtaVec.size(); ++i){
    float currentMass=massEtaVec.at(i);
    float currentEta=etaVec.at(i);
    bool found=false;
    float bestDiff=1000.;
    float matchedPhase=1000.;
    for(unsigned int j=0; j<massPhaseVec.size(); ++j){
      float currentPhaseMass=massPhaseVec.at(j);
      float currentDiff=fabs(currentMass-currentPhaseMass);
      if (currentDiff<matchAccuracy && currentDiff<bestDiff){
	bestDiff=currentDiff;
	matchedPhase=phaseVec.at(j);
	found=true;
      }
    }
    if (found){
      massResultVec.push_back(currentMass);
      etaResultVec.push_back(currentEta);
      phaseResultVec.push_back(matchedPhase);
    }
  }

  argandGraph= new TGraph();
  argandGraph->SetTitle("Argand Diagram #pi#pi #rightarrow #pi#pi");
  for(unsigned int i=0; i<massResultVec.size(); ++i){
    const std::complex<double> imagi(0.,1.);
    double etaVeci=etaResultVec.at(i);
    double phaseVeci=phaseResultVec.at(i)* TMath::Pi()/180.;
    complex<double> thePoint = 1./(2.*imagi)*(etaVeci*exp(2.*imagi*phaseVeci)-1.);
    argandGraph->SetPoint(i, thePoint.real(), thePoint.imag());
  }
  
  TExec *ex = new TExec("ex", "drawtext();");
  argandGraph->GetListOfFunctions()->Add(ex);
  argandGraph->Draw("ALP*");  
}

