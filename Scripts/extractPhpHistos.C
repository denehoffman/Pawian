#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TKey.h"
#include "TLegend.h"

#include <algorithm> 
#include <vector> 
//void extractPhpComplexHistos(std::string rootFileName);

void extractPhpHistos(std::string rootFileName="PhpFactor.root", double setyMin1d=-1.2, double setyMax1d=1){
 TFile* tFile= new TFile(rootFileName.c_str());
 
 TCanvas* canvasRealS = new TCanvas("RealS","RealS",20, 20, 900,900);
 canvasRealS->Divide(2,2);

 // TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
 // legend->AddEntry("","real","lpe");
 // legend->AddEntry("","imad","lpe");

 TH1F* phpAsnerImagHist = (TH1F*)tFile->Get("_phpAsnerImagHist");
 phpAsnerImagHist->Sumw2(false);
 phpAsnerImagHist->SetLineColor(kRed);
 phpAsnerImagHist->SetLineWidth(3);
 phpAsnerImagHist->SetStats(false);
 phpAsnerImagHist->SetMaximum(setyMax1d);
 phpAsnerImagHist->SetMinimum(setyMin1d);

 TH1F* phpAsnerRealHist = (TH1F*)tFile->Get("_phpAsnerRealHist");
 phpAsnerRealHist->Sumw2(false);
 phpAsnerRealHist->SetLineColor(kBlue);
 phpAsnerRealHist->SetLineWidth(3);
 phpAsnerRealHist->SetStats(false);
 phpAsnerRealHist->SetMaximum(setyMax1d);
 phpAsnerRealHist->SetMinimum(setyMin1d);

 TH1F* phpPenningtonImagHist = (TH1F*)tFile->Get("_phpPenningtonImagHist");
 phpPenningtonImagHist->Sumw2(false);
 phpPenningtonImagHist->SetLineColor(kRed);
 phpPenningtonImagHist->SetLineWidth(3);
 phpPenningtonImagHist->SetStats(false);
 phpPenningtonImagHist->SetMaximum(setyMax1d);
 phpPenningtonImagHist->SetMinimum(setyMin1d);

 TH1F* phpPenningtonRealHist = (TH1F*)tFile->Get("_phpPenningtonRealHist");
 phpPenningtonRealHist->Sumw2(false);
 phpPenningtonRealHist->SetLineColor(kBlue);
 phpPenningtonRealHist->SetLineWidth(3);
 phpPenningtonRealHist->SetStats(false);
 phpPenningtonRealHist->SetMaximum(setyMax1d);
 phpPenningtonRealHist->SetMinimum(setyMin1d);

 TH1F* phpDefaultComplImagHist = (TH1F*)tFile->Get("_phpDefaultComplImagHist");
 phpDefaultComplImagHist->Sumw2(false);
 phpDefaultComplImagHist->SetLineColor(kRed);
 phpDefaultComplImagHist->SetLineWidth(3);
 phpDefaultComplImagHist->SetMaximum(setyMax1d);
 phpDefaultComplImagHist->SetStats(false);
 phpDefaultComplImagHist->SetMinimum(setyMin1d);

 TH1F* phpDefaultComplRealHist = (TH1F*)tFile->Get("_phpDefaultComplRealHist");
 phpDefaultComplRealHist->Sumw2(false);
 phpDefaultComplRealHist->SetLineColor(kBlue);
 phpDefaultComplRealHist->SetLineWidth(3);
 phpDefaultComplRealHist->SetStats(false);
 phpDefaultComplRealHist->SetMaximum(setyMax1d);
 phpDefaultComplRealHist->SetMinimum(setyMin1d);

 TH1F* phpReidImagHist = (TH1F*)tFile->Get("_phpReidImagHist");
 phpReidImagHist->Sumw2(false);
 phpReidImagHist->SetLineColor(kRed);
 phpReidImagHist->SetLineWidth(3);
 phpReidImagHist->SetMaximum(setyMax1d);
 phpReidImagHist->SetMinimum(setyMin1d);
 phpReidImagHist->SetStats(false);

 TH1F* phpReidRealHist = (TH1F*)tFile->Get("_phpReidRealHist");
 phpReidRealHist->Sumw2(false);
 phpReidRealHist->SetLineColor(kBlue);
 phpReidRealHist->SetLineWidth(3);
 phpReidRealHist->SetStats(false);

 canvasRealS->cd(1);
 phpAsnerImagHist->Draw();
 phpAsnerRealHist->Draw("same");
 // legend->Draw();

 canvasRealS->cd(2);
 phpPenningtonImagHist->Draw();
 phpPenningtonRealHist->Draw("same");

 canvasRealS->cd(3);
 phpDefaultComplImagHist->Draw();
 phpDefaultComplRealHist->Draw("same");

 canvasRealS->cd(4);
 phpReidImagHist->Draw();
 phpReidRealHist->Draw("same");

 TCanvas* canvasRealM = new TCanvas("RealM","RealM",20, 20, 900,900);
 canvasRealM->Divide(2,2);

 TH1F* phpMAsnerImagHist = (TH1F*)tFile->Get("_phpMAsnerImagHist");
 phpMAsnerImagHist->Sumw2(false);
 phpMAsnerImagHist->SetLineColor(kRed);
 phpMAsnerImagHist->SetLineWidth(3);
 phpMAsnerImagHist->SetMaximum(setyMax1d);
 phpMAsnerImagHist->SetMinimum(setyMin1d);
 phpMAsnerImagHist->SetStats(false);

 TH1F* phpMAsnerRealHist = (TH1F*)tFile->Get("_phpMAsnerRealHist");
 phpMAsnerRealHist->Sumw2(false);
 phpMAsnerRealHist->SetLineColor(kBlue);
 phpMAsnerRealHist->SetLineWidth(3);
 phpMAsnerRealHist->SetStats(false);
 phpMAsnerRealHist->SetMaximum(setyMax1d);
 phpMAsnerRealHist->SetMinimum(setyMin1d);

 TH1F* phpMPenningtonImagHist = (TH1F*)tFile->Get("_phpMPenningtonImagHist");
 phpMPenningtonImagHist->Sumw2(false);
 phpMPenningtonImagHist->SetLineColor(kRed);
 phpMPenningtonImagHist->SetLineWidth(3);
 phpMPenningtonImagHist->SetMaximum(setyMax1d);
 phpMPenningtonImagHist->SetMinimum(setyMin1d);
 phpMPenningtonImagHist->SetStats(false);

 TH1F* phpMPenningtonRealHist = (TH1F*)tFile->Get("_phpMPenningtonRealHist");
 phpMPenningtonRealHist->Sumw2(false);
 phpMPenningtonRealHist->SetLineColor(kBlue);
 phpMPenningtonRealHist->SetLineWidth(3);
 phpMPenningtonRealHist->SetMaximum(setyMax1d);
 phpMPenningtonRealHist->SetMinimum(setyMin1d);
 phpMPenningtonRealHist->SetStats(false);

 TH1F* phpMDefaultComplImagHist = (TH1F*)tFile->Get("_phpMDefaultComplImagHist");
 phpMDefaultComplImagHist->Sumw2(false);
 phpMDefaultComplImagHist->SetLineColor(kRed);
 phpMDefaultComplImagHist->SetLineWidth(3);
 phpMDefaultComplImagHist->SetStats(false);
 phpMDefaultComplImagHist->SetMaximum(setyMax1d+5.);
 phpMDefaultComplImagHist->SetMinimum(setyMin1d);

 TH1F* phpMDefaultComplRealHist = (TH1F*)tFile->Get("_phpMDefaultComplRealHist");
 phpMDefaultComplRealHist->Sumw2(false);
 phpMDefaultComplRealHist->SetLineColor(kBlue);
 phpMDefaultComplRealHist->SetLineWidth(3);
 phpMDefaultComplRealHist->SetStats(false);
 phpMDefaultComplRealHist->SetMaximum(setyMax1d);
 phpMDefaultComplRealHist->SetMinimum(setyMin1d);
 
TH1F* phpMReidImagHist = (TH1F*)tFile->Get("_phpMReidImagHist");
 phpMReidImagHist->Sumw2(false);
 phpMReidImagHist->SetLineColor(kRed);
 phpMReidImagHist->SetLineWidth(3);
 phpMReidImagHist->SetMaximum(setyMax1d);
 phpMReidImagHist->SetMinimum(setyMin1d);
 phpMReidImagHist->SetStats(false);

 TH1F* phpMReidRealHist = (TH1F*)tFile->Get("_phpMReidRealHist");
 phpMReidRealHist->Sumw2(false);
 phpMReidRealHist->SetLineColor(kBlue);
 phpMReidRealHist->SetLineWidth(3);
 phpMReidRealHist->SetStats(false);
 phpMReidRealHist->SetMaximum(setyMax1d);
 phpMReidRealHist->SetMinimum(setyMin1d);

 canvasRealM->cd(1);
 phpMAsnerImagHist->Draw();
 phpMAsnerRealHist->Draw("same");

 canvasRealM->cd(2);
 phpMPenningtonImagHist->Draw();
 phpMPenningtonRealHist->Draw("same");

 canvasRealM->cd(3);
 phpMDefaultComplImagHist->Draw();
 phpMDefaultComplRealHist->Draw("same");

 canvasRealM->cd(4);
 phpMReidImagHist->Draw();
 phpMReidRealHist->Draw("same");






 TCanvas* canvasComplexS = new TCanvas("ComplexS","ComplexS",20, 20, 900,900);
  canvasComplexS->Divide(2,3);

 TH2F*	phpPenningtonRealHist2=(TH2F*)tFile->Get("_phpPenningtonRealHist2");
 phpPenningtonRealHist2->SetStats(false);
 TH2F*	phpPenningtonImagHist2=(TH2F*)tFile->Get("_phpPenningtonImagHist2");
 phpPenningtonImagHist2->SetStats(false);

 TH2F*	phpReidRealHist2=(TH2F*)tFile->Get("_phpReidRealHist2");
 phpReidRealHist2->SetStats(false);
 TH2F*	phpReidImagHist2=(TH2F*)tFile->Get("_phpReidImagHist2");
 phpReidImagHist2->SetStats(false);

 TH2F* phpDefaultComplRealHist2 =(TH2F*)tFile->Get("_phpDefaultComplRealHist2");
 phpDefaultComplRealHist2->SetStats(false);

 TH2F* phpDefaultComplRealAbove1Hist2 = (TH2F*) phpDefaultComplRealHist2->Clone();
 phpDefaultComplRealAbove1Hist2->GetXaxis()->SetRangeUser(0.2, 3.0); 

 TH2F* phpDefaultComplImagHist2 =(TH2F*)tFile->Get("_phpDefaultComplImagHist2");
 phpDefaultComplImagHist2->GetXaxis()->SetRangeUser(0.2, 3.0); 
 phpDefaultComplImagHist2->SetStats(false);

 canvasComplexS->cd(1);
 phpPenningtonRealHist2->Draw("colz");

 canvasComplexS->cd(2);
 phpPenningtonImagHist2->Draw("colz");

 canvasComplexS->cd(3);
 phpDefaultComplRealAbove1Hist2->Draw("colz");

 canvasComplexS->cd(4);
 phpDefaultComplImagHist2->Draw("colz");

 canvasComplexS->cd(5);
 phpReidRealHist2->Draw("colz");

 canvasComplexS->cd(6);
 phpReidImagHist2->Draw("colz");
}

