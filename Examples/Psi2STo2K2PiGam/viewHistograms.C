#include <vector>
#include <map>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <TFile>
#include <TH1F>
#include <TNtuple>
#include <TCanvas>

bool printToPDF = false;

void viewHistograms(TString fname="bin/gcc-4.1.2/debug/link-static/Psi2STo2K2PiGam.root", TString hypname="hypname", TString option){

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

  bool drawswitch = false;

  if (option.Contains("makepdf")) {printToPDF = true;}

  histVectData.push_back(invKKDataHist);
  histVectData.push_back(invKPiDataHist);
  histVectData.push_back(invPiPiDataHist);
  histVectData.push_back(invKPiPiViaK892DataHist);
  histVectData.push_back(cosPsiDataHist);
  histVectData.push_back(cosK892DataHist);
  histVectData.push_back(cosK1430DataHist);
  histVectData.push_back(cosK1430ViaK892DataHist);

  histVectMc.push_back(invKKFittedHist);
  histVectMc.push_back(invKPiFittedHist);
  histVectMc.push_back(invPiPiFittedHist);
  histVectMc.push_back(invKPiPiViaK892FittedHist);
  histVectMc.push_back(cosPsiFittedHist);
  histVectMc.push_back(cosK892FittedHist);
  histVectMc.push_back(cosK1430FittedHist);
  histVectMc.push_back(cosK1430ViaK892FittedHist);

  TCanvas* cmain = new TCanvas("cmain","cmain",1400,600);
  cmain->Divide(4,2);
  for(int i=0; i<histVectData.size(); i++) {
    cmain->cd(i+1);
    histVectData[i]->SetLineWidth(3);
    histVectData[i]->SetLineColor(2);
    histVectData[i]->Draw("E");
    histVectMc[i]->SetLineWidth(3);
    histVectMc[i]->Draw("same");
    } 
  
    cout << endl;

  histVectData2d.push_back(KPivsKPiDataHist);
  histVectData2d.push_back(KKvsPiPiDataHist);
  histVectData2d.push_back(KPiPivsPiPiDataHist);
  histVectData2d.push_back(KKPivsKPiDataHist);

  histVectMc2d.push_back(KPivsKPiFittedHist);
  histVectMc2d.push_back(KKvsPiPiFittedHist);
  histVectMc2d.push_back(KPiPivsPiPiFittedHist);
  histVectMc2d.push_back(KKPivsKPiFittedHist);

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


  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TCanvas* c2 = new TCanvas("c2","c2",1200,800);
  TCanvas* c3 = new TCanvas("c3","c3",1200,800);
  TCanvas* c4 = new TCanvas("c4","c4",1200,800);

  for(int i=0; i<4; i++) {
    if(i==0) c1->cd();
    if(i==1) c2->cd();
    if(i==2) c3->cd();
    if(i==3) c4->cd();
    histVectData[i]->SetLineWidth(4);
    histVectData[i]->SetLineColor(2);
    histVectData[i]->Draw("E");
    histVectMc[i]->SetLineWidth(4);
    histVectData[i]->GetXaxis()->SetTitle("invariante Masse / GeV/c^{2}");
    histVectData[i]->GetYaxis()->SetTitle("Ereignisse / 20 MeV/c^{2}");
    histVectMc[i]->Draw("same");
    TString histname = histVectData[i]->GetName();
    // cout << histname << endl;
    if(printToPDF) {
      histVectData[i]->SetTitle("");
      histVectMc[i]->SetTitle("");
      gStyle->SetOptStat(0);
      if(i==0) c1->Print("pdfplots/pwa_"+hypname+"_"+histname+".pdf");
      if(i==1) c2->Print("pdfplots/pwa_"+hypname+"_"+histname+".pdf");
      if(i==2) c3->Print("pdfplots/pwa_"+hypname+"_"+histname+".pdf");
      if(i==3) c4->Print("pdfplots/pwa_"+hypname+"_"+histname+".pdf");
    }
  } 

}

 


