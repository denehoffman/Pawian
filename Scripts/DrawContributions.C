//************************************************************************//
//                                                                        //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)                      //
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
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include <algorithm>    // std::max

TCanvas* contribCanvas = new TCanvas();
bool withresid=false;
bool withData=true;
int binning=1;
double x_min=0.;
double x_max=0.;

void setBinning(int theBinning){binning=theBinning;}

void setXminXmax(double thexMin, double thexMax){
    x_min=thexMin;
    x_max=thexMax;
  }

void DrawContributions(std::string rootFileNameData, std::string rootFileNames, std::string name, bool withLegend = false, std::string legendNames = "");

void DrawContributions(std::string rootFileNames, std::string name, bool withLegend = false,
		      std::string legendNames = "");

void DrawContrib2Hists(std::string rootFileName1, std::string rootFileName2, std::string name, bool withLegend = false, std::string legendNames = "");   


void DrawContributions(std::string rootFileNameData, std::string rootFileNames, std::string name, bool withLegend = false, 
                          std::string legendNames = "") {
 
  int colors[9] = {861, 902, 800, 825, 869, 877, 814, 808, 862}; 

  std::vector<TH1F*> histVec;
  std::stringstream rootFileNamesStrStr(rootFileNames);
  vector<std::string> rootFileNamesVec;
  std::string currentStr;
  while (rootFileNamesStrStr >> currentStr) rootFileNamesVec.push_back(currentStr);
  vector<std::string>::iterator it;
  for(it = rootFileNamesVec.begin(); it != rootFileNamesVec.end(); ++it) {
    std::cout << (*it) << std::endl;
  }

  TLegend *legend = new TLegend(0.19, 0.59, 0.39, 0.89);
  vector<std::string> legendNamesVec;
  if (withLegend) {
    std::stringstream legendNamesStrStr(legendNames);
    while (legendNamesStrStr >> currentStr) legendNamesVec.push_back(currentStr);
    for (it = legendNamesVec.begin(); it != legendNamesVec.end(); ++it){
      std::cout << (*it) << std::endl;
    }
    if (rootFileNamesVec.size() != legendNamesVec.size()) {
      std::cout << "rootFileNamesVec.size()= " << rootFileNamesVec.size() 
          << " != legendNamesVec.size(): " << legendNamesVec.size() << std::endl;
      return;
    }
  }


TPad * pad0 = new TPad("pad0","This is pad0",0.0,0.95,1.0,0.0);                    
  pad0->SetFillColor(0);
  TPad * pad1 = new TPad("pad1","This is pad1",0.0,1.,1.,0.85);                         
  pad1->SetFillColor(0);                                                                      

  if(withresid){
  pad0->Draw();
  pad1->Draw();
  pad0->cd();
  }

  double dataScale=1.;
  
  vector<std::string> histoNamesVec;

  std::stringstream histoNamesStrStr(name);
  while (histoNamesStrStr >> currentStr) histoNamesVec.push_back(currentStr);

  if(withData){

  TH1F* DataHist;
  TFile* dataTFile=new TFile(rootFileNameData.c_str());

  int count = 0;

  for(it=histoNamesVec.begin(); it!=histoNamesVec.end(); ++it){

      string dummy=it->substr(3);
      dummy="Data"+dummy;

      if(count ==0) DataHist = (TH1F*) dataTFile->Get(dummy.c_str());
          else{
              DataHist->Add((TH1F*)dataTFile->Get(dummy.c_str()));
          }
       count++;
  }
  
  DataHist->Rebin(binning);
  if(x_min!=0. && x_max!=0.) DataHist->GetXaxis()->SetRangeUser(x_min, x_max);
  DataHist->SetLineColor(kBlack);
  DataHist->SetMarkerColor(kBlack);
  DataHist->SetMarkerStyle(20);
  DataHist->SetLineWidth(2);
  //DataHist->Draw("E_SAME");
  legend->AddEntry(DataHist, "data", "l");
  DataHist->Draw("E");
  DataHist->SetMinimum(0.);
  dataScale=DataHist->Integral();
  }
  
  histVec.resize(rootFileNamesVec.size());
  
  TH1F* hFit; 

  for (unsigned int id = 0; id < rootFileNamesVec.size(); ++id){
    std::string currentFileName = rootFileNamesVec.at(id);
    TFile* currentTFile = new TFile(currentFileName.c_str());
    //TH1F* currentHist = (TH1F*)currentTFile->Get(name.c_str());
    TH1F* currentHist;

    int count =0;

    for(it=histoNamesVec.begin(); it!=histoNamesVec.end(); ++it){

        if(count ==0) currentHist = (TH1F*) currentTFile->Get(it->c_str());
        else{
            currentHist->Add((TH1F*)currentTFile->Get(it->c_str()));
        }
        count++;

    }
    currentHist->Rebin(binning);
    if(x_min!=0.&&x_max!=0.){
      currentHist->GetXaxis()->SetRangeUser(x_min, x_max);
    }

    //histVec.push_back(currentHist);
    
    if(id==0) dataScale=dataScale/currentHist->Integral();
    if(rootFileNamesVec.size()<9){
    currentHist->SetLineColor(colors[id]);
    }else{
        currentHist->SetLineColor(kBlack + id +1);
    }
    if (id==0){ 
        currentHist->SetLineColor(kRed); 
        currentHist->SetMarkerColor(kRed); 
    }
    currentHist->SetLineWidth(2);
    histVec[id] = currentHist;
    currentHist->Scale(dataScale);
    //if (0 ==id) 
    //    currentHist->Draw(); 
    //else 
    currentHist->Draw("E_HIST_SAME");
    currentHist->SetMaximum(1.05*(std::max(currentHist->GetMaximum(), histVec[0]->GetMaximum())));
    if (withLegend) 
      legend->AddEntry(currentHist, legendNamesVec.at(id).c_str(), "l");
  }

  if (withLegend) {
      legend->SetFillColor(0);
      legend->SetLineColor(kGray);
      legend->SetBorderSize(0);
      legend->SetShadowColor(0);
      legend->SetCornerRadius(0.05);
      legend->Draw("ARC");
      legend->SetFillColor(kWhite);

      //legend->SetFillColor(0);
      //legend->SetBorderSize(1);
      //legend->Draw(); 
  }

  if(withresid){
    TH1F* hresid= (TH1F*) histVec.at(0)->Clone();
    for(int i=0; i<histVec.at(0)->GetNbinsX();i++){ 
          double chisq=0.;
          if(histVec.at(0)->GetBinError(i)!=0)chisq=(histVec.at(0)->GetBinContent(i)-histVec.at(1)->GetBinContent(i))/sqrt(histVec.at(0)->GetBinError(i)*histVec.at(0)->GetBinError(i)+histVec.at(1)->GetBinError(i)*histVec.at(1)->GetBinError(i));
          if(histVec.at(0)->GetBinError(i)==0) cout << i << endl;
          if(histVec.at(0)->GetBinContent(i)-histVec.at(1)->GetBinContent(i)==0) cout << i << endl;
          hresid->SetBinContent(i,chisq);
	  hresid->SetBinError(i,0.);

	  pad1->cd();
    
    hresid->SetMarkerStyle(20);
    hresid->SetMarkerSize(1.);
    TLine *line1 = new TLine(hresid->GetBinLowEdge(1),0,hresid->GetBinLowEdge(hresid->GetNbinsX()+1),0);
    TLine *line2 = new TLine(hresid->GetBinLowEdge(1),-3,hresid->GetBinLowEdge(hresid->GetNbinsX()+1),-3);
    TLine *line3 = new TLine(hresid->GetBinLowEdge(1),3.,hresid->GetBinLowEdge(hresid->GetNbinsX()+1),3.);
    
    line1->SetLineColor(kRed);  line1->SetLineWidth(2);
    line2->SetLineColor(kRed);  line2->SetLineWidth(2);
    line3->SetLineColor(kRed);  line3->SetLineWidth(2);
    line1->SetLineStyle(7);
        
    hresid->SetTitle("");
    hresid->Draw("Hist P");
    hresid->GetYaxis()->SetRangeUser(-3.5,3.5);
    hresid->GetYaxis()->SetTitle("Pull");
    hresid->GetYaxis()->SetLabelSize(0.1); 
    hresid->GetYaxis()->SetTitleSize(0.15);
    hresid->GetYaxis()->SetTitleOffset(0.25);
    line1->Draw();
    line2->Draw();
    line3->Draw();
    }
  }
}

void DrawContributions(std::string rootFileNames, std::string name, bool withLegend,
		       std::string legendNames){
  std::stringstream rootFileNamesStrStr(rootFileNames);
  std::string firstRootFileName;
  rootFileNamesStrStr >> firstRootFileName;
  DrawContributions(firstRootFileName, rootFileNames, name, withLegend,legendNames);
}


  void DrawContrib2Hists(std::string rootFileName1, std::string rootFileName2, std::string name, bool withLegend, std::string legendNames){
  std::string rootFileNames=rootFileName1+" "+rootFileName2;
  withresid=true;
  withData=false;
  DrawContributions(rootFileNames, name, withLegend, legendNames);
}
