#include <vector>
#include <map>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <TFile>
#include <TH1F>
#include <TNtuple>


void viewHistograms(TString fname="bin/gcc-4.1.2/debug/link-static/Psi2STo2K2PiGam.root"){

  using namespace std;
  gROOT->SetStyle("Plain");

  gStyle->SetCanvasColor(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111); 

  TFile* f1 = new TFile(fname,"READ");

  std::vector<TH1F*> histVectData;
  std::vector<TH1F*> histVectMc;
  std::vector<TH2F*> histVectData2d;
  std::vector<TH2F*> histVectMc2d;

  histVectData.push_back(_invKKDataHist);
  histVectData.push_back(_invKPiDataHist);
  histVectData.push_back(_invPiPiDataHist);
  histVectData.push_back(_invKPiPiViaK892DataHist);
  histVectData.push_back(_cosPsiDataHist);
  histVectData.push_back(_cosK892DataHist);
  histVectData.push_back(_cosK1430DataHist);
  histVectData.push_back(_cosK1430ViaK892DataHist);

  histVectMc.push_back(_invKKFittedHist);
  histVectMc.push_back(_invKPiFittedHist);
  histVectMc.push_back(_invPiPiFittedHist);
  histVectMc.push_back(_invKPiPiViaK892FittedHist);
  histVectMc.push_back(_cosPsiFittedHist);
  histVectMc.push_back(_cosK892FittedHist);
  histVectMc.push_back(_cosK1430FittedHist);
  histVectMc.push_back(_cosK1430ViaK892FittedHist);

  TCanvas* cmain = new TCanvas("cmain","cmain",1400,600);
  cmain->Divide(4,2);
  for(int i=0; i<histVectData.size(); i++) {
    cmain->cd(i+1);
    histVectData[i]->SetLineColor(2);
    histVectData[i]->Draw("E");
    histVectMc[i]->Draw("same");
    } 
  
    cout << endl;


  histVectData2d.push_back(_KPivsKPiDataHist2d);
  histVectData2d.push_back(_KKvsPiPiDataHist2d);
  histVectData2d.push_back(_KPiPivsPiPiDataHist2d);
  histVectData2d.push_back(_KKPivsKPiDataHist2d);

  histVectMc2d.push_back(_KPivsKPiFittedHist2d);
  histVectMc2d.push_back(_KKvsPiPiFittedHist2d);
  histVectMc2d.push_back(_KPiPivsPiPiFittedHist2d);
  histVectMc2d.push_back(_KKPivsKPiFittedHist2d);

  TCanvas* cmain2d = new TCanvas("cmain2d","cmain2d",1400,600);
  cmain2d->Divide(4,2);
  for(int i=0; i<histVectData2d.size(); i++) {
    cmain2d->cd(i+1);
    histVectData2d[i]->Draw("colz");
    cmain2d->cd(i+1+histVectData2d.size());
    histVectMc2d[i]->SetMaximum(histVectData2d[i]->GetMaximum());
    histVectMc2d[i]->Draw("colz");
    } 
  
    cout << endl;




}

 


