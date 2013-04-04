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
#include "Examples/Tutorial/MinuitFit/MinuitFit.hh"
#include "Examples/Tutorial/MinuitFit/MinuitFitFcn.hh"

#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include "ErrLogger/ErrLogger.hh"

MinuitFit::MinuitFit(double p0, double p1, double p2, double p3, double sigma) :
  _theTFile(0),
  _sigma(sigma)
{

  // Display parameters for test distribution 

  cout << endl;
  Info <<"Set p0 as "<< p0 << endmsg;
  Info <<"Set p1 as "<< p1 << endmsg;
  Info <<"Set p2 as "<< p2 << endmsg;
  Info <<"Set p3 as "<< p3 << endmsg;
  Info <<"Set sigma as "<< sigma << endmsg;

  // Generate test distribution and smear them with a gaussian
  TRandom randomNumber;
  for(int i=0; i<1000; i++) {
    double tmpXvalue=static_cast<double>((rand() % 10000 + 1)) / 100;
    double tmpYvalue=randomNumber.Gaus(p0 + p1 * tmpXvalue + p2 * tmpXvalue * tmpXvalue + p3 * tmpXvalue * tmpXvalue * tmpXvalue, _sigma);
    _xValue.push_back(tmpXvalue);
    _yValue.push_back(tmpYvalue);
  }

}

double MinuitFit::calcChiSqr(const std::vector<double>& minPar){
 
  // Calculate chi^2 for current set of fit parameters 
  double result=0.;
  for (unsigned int i=0; i<_xValue.size(); i++){
    double yValFit=minPar[0]+minPar[1]*_xValue[i]+minPar[2]*_xValue[i]*_xValue[i]+minPar[3]*_xValue[i]*_xValue[i]*_xValue[i];
    double yValExp=_yValue[i];
    double tmpChi=((yValExp-yValFit)*(yValExp-yValFit))/(_sigma*_sigma);
    result+=tmpChi;
  }
  return result;
}


void MinuitFit::drawGraph(double a, double b, double c, double d){

  // Create arrays (which will be needed for feeding the TGraph)
  const unsigned int dataEntries=_xValue.size();
  double x[dataEntries];
  double y[dataEntries]; 
  double xerr[dataEntries];
  double yerr[dataEntries]; 

  // Convert vectors to arrays
  for(unsigned int i=0; i<dataEntries; i++) {
    x[i]=_xValue[i];
    y[i]=_yValue[i];
    xerr[i] = 0;
    yerr[i] = _sigma;
  }

  // Create ROOT file                                                                                                             
  _theTFile = new TFile("myFit.root","RECREATE");

  // Create TGraph and write to ROOT file 
  TCanvas* c1=new TCanvas("c1","c1",1200,800);
  c1->cd();
  TGraphErrors* linDist = new TGraphErrors(100, x, y, xerr, yerr);
  linDist->SetMarkerStyle(6);
  linDist->Draw("AP");
  linDist->Write("myGraph");

  // Draw fit with calculated parameters
  TF1* fit1 = new TF1("fit1","pol3",0.,100.);
  fit1->SetParameters(a, b, c, d);
  fit1->SetLineColor(46);
  fit1->SetLineWidth(3);
  fit1->Draw();
  fit1->Write();

}


MinuitFit::~MinuitFit()
{
   _theTFile->Write();
   _theTFile->Close();
}

