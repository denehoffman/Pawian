#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TROOT.h"
#include "TStyle.h"
#include <Riostream.h>
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

  TLorentzVector* KsVec=new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector* pionpVec=new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector* pionmVec=new TLorentzVector(0.,0.,0.,0.);

class Histos{
public:
  Histos(std::string suffix){
    std::string currentName="K_{S} #pi^{+} #pi^{-} Dpl"+suffix;
    kSPipPimDplH2=new TH2F("kSPipPimDplH2",currentName.c_str(), 75, 0.45*0.45+0.14*0.14, 3.5, 75, 0.45*0.45+0.14*0.14, 3.5);

    currentName="M(#pi^{+} #pi^{-}) H1"+suffix;
    invPiPiH1=new TH1F("invPiPiH1",currentName.c_str(), 200, 0.25, 1.5);
    invPiPiH1->Sumw2();
    invPiPiH1->SetLineWidth(2);

    currentName="M(K_{S} #pi^{+}) H1"+suffix;
    invKsPipH1=new TH1F("invKsPipH1",currentName.c_str(), 200, 0.6, 1.9);
    invKsPipH1->Sumw2();
    invKsPipH1->SetLineWidth(2);

    currentName="M(K_{S} #pi^{-}) H1"+suffix;
    invKsPimH1=new TH1F("invKsPimH1",currentName.c_str(), 200, 0.6, 1.9);
    invKsPimH1->Sumw2();
    invKsPimH1->SetLineWidth(2);
  }
  ~Histos(){;}
  void scale(double factor){
    kSPipPimDplH2->Scale(factor);
    invPiPiH1->Scale(factor);
    invKsPipH1->Scale(factor);
    invKsPimH1->Scale(factor);
  }
  TH2F* kSPipPimDplH2;
  TH1F* invPiPiH1;
  TH1F* invKsPipH1;
  TH1F* invKsPimH1;
};

Float_t plotDpl( TTree* theTree, Histos &theHistos){
  Float_t weight=0;
  Float_t noOfWeightedEvts=0.;

  theTree->SetBranchAddress("Kshort", &KsVec);
  theTree->SetBranchAddress("pion+", &pionpVec);
  theTree->SetBranchAddress("pion-", &pionmVec);
  theTree->SetBranchAddress("weight", &weight);

  Int_t noOfEntries=theTree->GetEntries();
  for (Int_t i=0; i<noOfEntries; ++i){
    theTree->GetEntry(i);
    TLorentzVector currentKsPip= *KsVec+*pionpVec;
    TLorentzVector currentKsPim= *KsVec+*pionmVec;
    TLorentzVector currentPiPi= *pionmVec+*pionpVec;


    theHistos.kSPipPimDplH2->Fill(currentKsPip.M()*currentKsPip.M(), currentKsPim.M()*currentKsPim.M(), weight);
    theHistos.invPiPiH1->Fill(currentPiPi.M(), weight);
    theHistos.invKsPipH1->Fill(currentKsPip.M(), weight);
    theHistos.invKsPimH1->Fill(currentKsPim.M(), weight);
    
    noOfWeightedEvts+=weight;
    
  }
  return noOfWeightedEvts;
  //scale histo
  //  theHist->Scale(10000./noOfEntries);
}

void drawResult(){
  TFile* f=new TFile("pawianHists_res.root");
  f->ls();
  TTree*  dataTree=(TTree*) f->Get("_dataFourvecs");
  Histos histosData("Data");
  Histos histosFit("Fit");

  dataTree->Print();

  TTree*  fitTree=(TTree*) f->Get("_fittedFourvecs");
  fitTree->Print();

  Float_t noOfWeightedDataEvts=plotDpl(dataTree, histosData);
  Float_t noOfWeightedFitEvts=plotDpl(fitTree, histosFit);

  histosFit.scale(noOfWeightedDataEvts/noOfWeightedFitEvts);
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c1=new TCanvas();
  //  c1->Size(200.,200.);
  c1->SetFillColor(0);
  c1->Divide(2,2);
  c1->cd(1);
  histosData.kSPipPimDplH2->Draw("colz");
  c1->cd(2);
  histosFit.kSPipPimDplH2->Draw("colz");
  

  TCanvas* c2=new TCanvas();
  //  c2->Size(200.,200.);
  c2->SetFillColor(0);
  c2->Divide(2,2);
  c2->cd(1);
  histosData.invPiPiH1->Draw();
  histosFit.invPiPiH1->SetLineColor(kRed);
  histosFit.invPiPiH1->Draw("same");
  c2->cd(2);
  histosData.invKsPipH1->Draw();
  histosFit.invKsPipH1->SetLineColor(kRed);
  histosFit.invKsPipH1->Draw("same");
  c2->cd(3);
  histosData.invKsPimH1->Draw();
  histosFit.invKsPimH1->SetLineColor(kRed);
  histosFit.invKsPimH1->Draw("same");
}
