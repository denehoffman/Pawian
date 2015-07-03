#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include <algorithm>    // std::max

void DrawContibutions(std::string rootFileNames, std::string name, bool withLegend=false, std::string legendNames="");

void DrawContibutions(std::string rootFileNames, std::string name, bool withLegend, std::string legendNames){
  //  std::vector<TH1F*> histVec;
  TH1F* histVec[20];
  
  std::stringstream rootFileNamesStrStr(rootFileNames);
  vector<std::string> rootFileNamesVec;
  std::string currentStr;
  while (rootFileNamesStrStr >> currentStr) rootFileNamesVec.push_back(currentStr);
  vector<std::string>::iterator it;
  for(it=rootFileNamesVec.begin(); it!=rootFileNamesVec.end(); ++it){
    std::cout << (*it) << std::endl;
  }

  TLegend *legend = new TLegend(0.19,0.59,0.39,0.89);
  vector<std::string> legendNamesVec;
  if(withLegend){
    std::stringstream legendNamesStrStr(legendNames);
    while (legendNamesStrStr >> currentStr) legendNamesVec.push_back(currentStr);
    for(it=legendNamesVec.begin(); it!=legendNamesVec.end(); ++it){
      std::cout << (*it) << std::endl;
    }
    if(rootFileNamesVec.size() != legendNamesVec.size()){
      std::cout << "rootFileNamesVec.size()= " << rootFileNamesVec.size() <<" != legendNamesVec.size(): " << legendNamesVec.size() << std::endl;
      return;
    }
  }

  for (unsigned int id=0; id<rootFileNamesVec.size(); ++id){
    std::string currentFileName=rootFileNamesVec.at(id);
    TFile* currentTFile=new TFile(currentFileName.c_str());
    TH1F* currentHist=(TH1F*) currentTFile->Get(name.c_str());
    currentHist->SetLineColor(kBlack+id);
    currentHist->SetLineWidth(2);
    histVec[id]=currentHist;
    if (id==0) currentHist->Draw(); 
    else currentHist->Draw("same");
    currentHist->SetMaximum(std::max(currentHist->GetMaximum(), histVec[0]->GetMaximum()));
    if (withLegend) legend->AddEntry(currentHist, legendNamesVec.at(id).c_str(),"l");
  }

  if(withLegend){
    legend->SetFillColor(0);
    legend->SetBorderSize(1);
    legend->Draw(); 
  } 
}
