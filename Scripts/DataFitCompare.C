#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TKey.h"
#include <algorithm>    // std::max

void DataFitCompare(std::string rootFileName, std::string histname, bool createcanvas=true);
double CalcHistoChi2(TH1F* histoData, TH1F* histoFit);
void DataFitComp2(TH1F* histoData, TH1F* histoFit);
void DataFitCompareAll(std::string rootFileName, bool saveImage=false);

void DataFitComp2(TH1F* histoData, TH1F* histoFit){


   histoData->SetLineColor(kRed);
   histoFit->SetLineColor(kBlack);

   histoData->SetLineWidth(2);
   histoFit->SetLineWidth(2);

   histoFit->SetFillColor(kYellow-10);

   histoFit->SetMinimum(0);
   histoFit->SetMaximum(histoFit->GetMaximum()*1.3);

   histoFit->Draw("hist");
   histoFit->Draw("same");
   histoData->Draw("same");
}

void DataFitCompare(std::string rootFileName, std::string histname, bool createcanvas){
  TFile* tFile= new TFile(rootFileName.c_str());
   std::string dataName = "Data" + histname;
   std::string fitName = "Fit" + histname;

   TH1F* histoData = (TH1F*)tFile->Get(dataName.c_str());
   TH1F* histoFit = (TH1F*)tFile->Get(fitName.c_str());

   if(createcanvas)
     TCanvas* canvas = new TCanvas("canvas","c1",1000,1000);

   CalcHistoChi2(histoData, histoFit);

   DataFitComp2(histoData, histoFit);
}

double CalcHistoChi2(TH1F* histoData, TH1F* histoFit){

   int nbins1 = histoData->GetNbinsX();
   int nbins2 = histoFit->GetNbinsX();

   if(nbins1 != nbins2){
      std::cout << "ERROR: nbins1!=nbins2" << std::endl;
   }

   double chi2=0;
   int nonzerobins=0;
   for(int i=1;i<=nbins1;i++){
      if(histoData->GetBinContent(i) == 0 || histoFit->GetBinContent(i) == 0)
	 continue;

      double error = sqrt(histoData->GetBinContent(i) + histoFit->GetBinContent(i));

      chi2+=pow(histoData->GetBinContent(i) - histoFit->GetBinContent(i), 2)/(error*error);
      //      std::cout << pow(histoData->GetBinContent(i) - histoFit->GetBinContent(i), 2)/(error*error) << std::endl;
      nonzerobins++;
   }
   
   std::cout << "Chi2 for histo " << histoData->GetName() << "\t" << chi2/nonzerobins 
	     << "p =" << TMath::Prob(chi2, nonzerobins) << std::endl;


   return chi2;
}

void DataFitCompareAll(std::string rootFileName, bool saveImage){
  TFile* tFile= new TFile(rootFileName.c_str());   
   std::vector<std::string> histnames;

   TIter nextkey(tFile->GetListOfKeys());
   TKey *key;

   while (key = (TKey*)nextkey()) {
      
      TH1F *histo = (TH1F*)key->ReadObj();
      const char* name = histo->GetName();
      //      std::cout << name << std::endl;
      std::string sname(name);
      std::string pname = sname.substr(0, 4);
      std::string aname = sname.substr(4);
      
      if(pname == "Data"){
	 histnames.push_back(aname);
      }
   }

   TCanvas* canvas = new TCanvas("canvas","c1",1000,1000);
   int numhists = histnames.size();

   canvas->Divide(std::ceil(sqrt(numhists)), std::ceil(sqrt(numhists)));

   for(int i=0; i<numhists; i++){
      canvas->cd(i+1);
      DataFitCompare(rootFileName, histnames.at(i), false);
   }

   if(saveImage){
      canvas->SaveAs("DataFitCompFile.png");
   }
}

