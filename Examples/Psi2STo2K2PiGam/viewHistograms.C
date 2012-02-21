#include <vector>
#include <map>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <TFile>
#include <TH1F.h>
#include <TH2F.h>
#include <TNtuple>
#include <TCanvas>
#include <TF1.h>

TNtuple* ntdata;

float datainvmassK1pi1;
float datainvmassK1pi2;
float datainvmassK2pi1;
float datainvmassK2pi2;
float datainvmassK1pi1pi2;
float datainvmassK2pi1pi2;
float datainvmasspi1pi2;
float datainvmassK1K2pi1;
float datainvmassK1K2pi2;
float datainvmassK1K2;

float datacosthetapi1;
float datacosthetapi2;
float datacosthetapipiViaK1pi1pi2;
float datacosthetapipiViaK2pi1pi2;
float datacosthetaK1K2ViaK1K2pi1;
float datacosthetaK1K2ViaK1K2pi2;

std::vector<float> vecdatainvmassK1pi1;
std::vector<float> vecdatainvmassK1pi2;
std::vector<float> vecdatainvmassK2pi1;
std::vector<float> vecdatainvmassK2pi2;
std::vector<float> vecdatacosthetapi1;
std::vector<float> vecdatacosthetapi2;

std::vector<float> vecdatainvmassK1pi1pi2;
std::vector<float> vecdatainvmassK2pi1pi2;
std::vector<float> vecdatainvmasspi1pi2;
std::vector<float> vecdatacosthetapipiViaK1pi1pi2;
std::vector<float> vecdatacosthetapipiViaK2pi1pi2;

std::vector<float> vecdatainvmassK1K2pi1;
std::vector<float> vecdatainvmassK1K2pi2;
std::vector<float> vecdatainvmassK1K2;
std::vector<float> vecdatacosthetaK1K2ViaK1K2pi1;
std::vector<float> vecdatacosthetaK1K2ViaK1K2pi2;

float datacosThetapi1;
float datacosThetapi2;

Int_t numOfEntriesData;

TH2F* datainvMassKpivsKpi;
TH1F* datainvMassKpiangular;

TH2F* datainvMassKpipivspipi;
TH1F* datainvMassKpipiangular;

bool printToPDF = false;

void viewHistograms(TString fname="bin/gcc-4.1.2/debug/link-static/Psi2STo2K2PiGam.root", TString hypname="hypname", TString option){

  using namespace std;
  gROOT->SetStyle("Plain");
  
  gStyle->SetCanvasColor(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111); 
  
  TFile* f1 = new TFile(fname,"READ");
  ntdata = (TNtuple*)f1->Get("dataTuple");
  numOfEntriesData = ntdata->GetEntries();

  ntdata->SetBranchAddress("mk1pi1", &datainvmassK1pi1);
  ntdata->SetBranchAddress("mk1pi2", &datainvmassK1pi2);
  ntdata->SetBranchAddress("mk2pi1", &datainvmassK2pi1);
  ntdata->SetBranchAddress("mk2pi2", &datainvmassK2pi2);
  ntdata->SetBranchAddress("costhetapi1", &datacosthetapi1);
  ntdata->SetBranchAddress("costhetapi2", &datacosthetapi2);

  ntdata->SetBranchAddress("mk1pi1pi2", &datainvmassK1pi1pi2);
  ntdata->SetBranchAddress("mk2pi1pi2", &datainvmassK2pi1pi2);
  ntdata->SetBranchAddress("mpipi", &datainvmasspi1pi2);
  ntdata->SetBranchAddress("costhetapipiViaK1pipi", &datacosthetapipiViaK1pi1pi2);
  ntdata->SetBranchAddress("costhetapipiViaK2pipi", &datacosthetapipiViaK2pi1pi2);

  for(int i = 0; i<numOfEntriesData; i++) {
    ntdata->GetEntry(i);
    vecdatainvmassK1pi1.push_back(datainvmassK1pi1);
    vecdatainvmassK1pi2.push_back(datainvmassK1pi2);
    vecdatainvmassK2pi1.push_back(datainvmassK2pi1);
    vecdatainvmassK2pi2.push_back(datainvmassK2pi2);
    vecdatacosthetapi1.push_back(datacosthetapi1);
    vecdatacosthetapi2.push_back(datacosthetapi2);

    vecdatainvmassK1pi1pi2.push_back(datainvmassK1pi1pi2);
    vecdatainvmassK2pi1pi2.push_back(datainvmassK2pi1pi2);
    vecdatainvmasspi1pi2.push_back(datainvmasspi1pi2);
    vecdatacosthetapipiViaK1pi1pi2.push_back(datacosthetapipiViaK1pi1pi2);
    vecdatacosthetapipiViaK2pi1pi2.push_back(datacosthetapipiViaK2pi1pi2);
  }

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

  datainvMassKpivsKpi = new TH2F("datainvMassKpivsKpi", "datainvMassKpivsKpi", 48, 0.5, 2.9, 48, 0.5, 2.9);
  datainvMassKpiAngular = new TH1F("datainvMassKpiAngular", "datainvMassKpiAngular", 20, -1., 1.);

  datainvMassKpipivspipi = new TH2F("datainvMassKpipivspipi", "datainvMassKpipivspipi", 46, 0.2, 2.5, 48, 0.8, 3.2);
  datainvMassKpipiAngular = new TH1F("datainvMassKpipiAngular", "datainvMassKpipiAngular", 20, -1., 1.);

  for(i=0; i<numOfEntriesData; i++) {
    if( (vecdatainvmassK1pi1[i]>0.870 && vecdatainvmassK1pi1[i]<0.930) && (vecdatainvmassK2pi2[i]>1.8 && vecdatainvmassK2pi2[i]<2.3) )
      {
	datainvMassKpivsKpi->Fill(vecdatainvmassK1pi1[i],vecdatainvmassK2pi2[i]);
	datainvMassKpivsKpi->Fill(vecdatainvmassK2pi2[i],vecdatainvmassK1pi1[i]);
	datainvMassKpiAngular->Fill(vecdatacosthetapi2[i]);
      }
    if( (vecdatainvmassK1pi2[i]>0.870 && vecdatainvmassK1pi2[i]<0.930) && (vecdatainvmassK2pi1[i]>1.8 && vecdatainvmassK2pi1[i]<2.3) )
      {
	datainvMassKpivsKpi->Fill(vecdatainvmassK1pi2[i],vecdatainvmassK2pi1[i]);
	datainvMassKpivsKpi->Fill(vecdatainvmassK2pi1[i],vecdatainvmassK1pi2[i]);
	datainvMassKpiAngular->Fill(vecdatacosthetapi1[i]);
      }
    if( (vecdatainvmassK2pi1[i]>0.870 && vecdatainvmassK2pi1[i]<0.930) && (vecdatainvmassK1pi2[i]>1.8 && vecdatainvmassK1pi2[i]<2.3) )
      {
	datainvMassKpivsKpi->Fill(vecdatainvmassK2pi1[i],vecdatainvmassK1pi2[i]);
	datainvMassKpivsKpi->Fill(vecdatainvmassK1pi2[i],vecdatainvmassK2pi1[i]);
	datainvMassKpiAngular->Fill(vecdatacosthetapi2[i]);
      }
    if( (vecdatainvmassK2pi2[i]>0.870 && vecdatainvmassK2pi2[i]<0.930) && (vecdatainvmassK1pi1[i]>1.8 && vecdatainvmassK1pi1[i]<2.3) )
      {
	datainvMassKpivsKpi->Fill(vecdatainvmassK2pi2[i],vecdatainvmassK1pi1[i]);
	datainvMassKpivsKpi->Fill(vecdatainvmassK1pi1[i],vecdatainvmassK2pi2[i]);
	datainvMassKpiAngular->Fill(vecdatacosthetapi1[i]);
      }
  }       

  TCanvas* cangular1 = new TCanvas("cangular1","cangular1",1200,600);
  cangular1->Divide(2,1);
  cangular1->cd(1);
  datainvMassKpivsKpi->Draw("colz");
  cangular1->cd(2);
  datainvMassKpiAngular->Draw();
  
 //  TF1* fit1 = new TF1("fit1","gaus(0)",-10.,10.);
//   fit1->SetParameters(5,histofpull->GetMean(),histofpull->GetRMS());
//   fit1->SetLineColor(4);
//   fit1->SetLineWidth(3);
//   datainvMassKpiAngular->Fit("fit1","LMRS","",-10.,10.);

  float lower_limit = 2.3; // 2.3;
  float upper_limit = 2.6; // 2.6;
    
  for(i=0; i<numOfEntriesData; i++) {
    if( (vecdatainvmasspi1pi2[i] > 0.94 && vecdatainvmasspi1pi2[i] < 1.01 ) && (vecdatainvmassK1pi1pi2[i] > lower_limit && vecdatainvmassK1pi1pi2[i] < upper_limit) )
      {
	datainvMassKpipivspipi->Fill(vecdatainvmasspi1pi2[i], vecdatainvmassK1pi1pi2[i]);
	datainvMassKpipiAngular->Fill(vecdatacosthetapipiViaK1pi1pi2[i]);
      }
    if( (vecdatainvmasspi1pi2[i] > 0.94 && vecdatainvmasspi1pi2[i] < 1.01 ) && (vecdatainvmassK2pi1pi2[i] > lower_limit && vecdatainvmassK2pi1pi2[i] < upper_limit) )
      {
	datainvMassKpipivspipi->Fill(vecdatainvmasspi1pi2[i], vecdatainvmassK2pi1pi2[i]);
	datainvMassKpipiAngular->Fill(vecdatacosthetapipiViaK2pi1pi2[i]);
      }

  } 

  TCanvas* cangular2 = new TCanvas("cangular2","cangular2",1200,600);
  cangular2->Divide(2,1);
  cangular2->cd(1);
  datainvMassKpipivspipi->Draw("colz");
  cangular2->cd(2);
  datainvMassKpipiAngular->Draw();

  TF1* fit2 = new TF1("fit2","[0]+pow([1]*(3/2*x*x-1/2)+[2],2)",-1.,1.);
  fit2->SetParameters(3,5,0.3);
  fit2->SetLineColor(4);
  fit2->SetLineWidth(3);
  datainvMassKpipiAngular->Fit("fit2","LMRS","",-1.,1.);

}

 


