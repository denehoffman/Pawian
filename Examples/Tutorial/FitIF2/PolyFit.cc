#include <getopt.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include "Examples/Tutorial/FitIF2/PolyFit.hh"

#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include "ErrLogger/ErrLogger.hh"

PolyFit::PolyFit(double p0, double p1, double p2, double p3, double sigma) :
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

PolyFit::PolyFit(const PolyFit& inFit) :
  _theTFile(0)
{

  _myGraph = inFit.getGraph();

  _xValue = inFit.getX();
  _yValue = inFit.getY();

  _sigma = inFit.getSig();
}

double PolyFit::evaluator(){
 
  // Calculate chi^2 for current set of fit parameters 
  double result=0.;
  for (unsigned int i=0; i<_xValue.size(); i++){
    double yValFit=parValue_[0]+parValue_[1]*_xValue[i]+parValue_[2]*_xValue[i]*_xValue[i]+parValue_[3]*_xValue[i]*_xValue[i]*_xValue[i];
    double yValExp=_yValue[i];
    double tmpChi=((yValExp-yValFit)*(yValExp-yValFit))/(_sigma*_sigma);
    result+=tmpChi;
  }
  return result;
}

void PolyFit::setStartPar(int num, double *par, double* min, double* max, double* err){
 
  for(int i=0; i<num ; i++){
    parValue_.push_back(par[i]);
    minValue_.push_back(min[i]);
    maxValue_.push_back(max[i]);
    parError_.push_back(err[i]);
  }

}

void PolyFit::drawGraph(double a, double b, double c, double d){

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


PolyFit::~PolyFit()
{
   _theTFile->Write();
   _theTFile->Close();
}

