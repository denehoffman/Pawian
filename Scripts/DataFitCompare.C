//************************************************************************//
//									  //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TKey.h"
#include <algorithm>    // std::max

double CalcHistoChi2(TH1F* histoData, TH1F* histoFit);
double CalcKolmogorov(TH1F* histoData, TH1F* histoFit);
void DataFitCompare(std::string rootFileName, std::string histname, bool createcanvas=true);
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
     //     TCanvas* canvas = new TCanvas("canvas","c1",1000,1000);

     std::cout << "\n" << std::endl;
   CalcHistoChi2(histoData, histoFit);
   CalcKolmogorov(histoData, histoFit);
   std::cout << "\n" << std::endl;
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
      // if(histoData->GetBinContent(i) == 0 && histoFit->GetBinContent(i) == 0)
      // 	 continue;


      //     double error = sqrt(histoData->GetBinContent(i) + histoFit->GetBinContent(i));
      double error =  sqrt(histoData->GetBinError(i)*histoData->GetBinError(i) + histoFit->GetBinError(i)*histoFit->GetBinError(i));

      chi2+=pow(histoData->GetBinContent(i) - histoFit->GetBinContent(i), 2)/(error*error);
      //      std::cout << pow(histoData->GetBinContent(i) - histoFit->GetBinContent(i), 2)/(error*error) << std::endl;
      nonzerobins++;
   }
   double chi2Root=histoData->Chi2Test(histoFit,"WW");
    std::cout << "Chi2 for histo " << histoData->GetName() << "\t" << chi2/nonzerobins 
	     << "\tp =" << TMath::Prob(chi2, nonzerobins) 
	      << "\tpRoot=" << chi2Root  << std::endl;


   return chi2;
}

double CalcKolmogorov(TH1F* histoData, TH1F* histoFit){
  double ksTestResult = histoData->KolmogorovTest(histoFit);
  std::cout << "Kolmogorov-Smirnov test for histo " << histoData->GetName() 
	    << "\tp = " << ksTestResult << std::endl;
  return ksTestResult;  
}

void DataFitCompareAll(std::string rootFileName, bool saveImage){
  TFile* tFile= new TFile(rootFileName.c_str());   
   std::vector<std::string> histnamesMasses;
   std::vector<std::string> histnamesAnglesHeli;
   std::vector<std::string> histnamesAnglesGJ;

   std::string strTheta("Theta");
   std::string strPhi("Phi");
   std::string strGJ("GJ");

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
	//sort angle hists
	if ( (aname.find("Theta") != string::npos) || (aname.find("Phi") != string::npos) ){
	  //angle hists in GJ system
	  if ( aname.find("GJ") != string::npos) histnamesAnglesGJ.push_back(aname);
	  else histnamesAnglesHeli.push_back(aname);      
	}
	else histnamesMasses.push_back(aname);
      }
   }

   TCanvas* canvasMass = new TCanvas("Masses","Masses",20, 20, 900,900);
   int numhistsMasses = histnamesMasses.size();
   canvasMass->Divide(std::ceil(sqrt(numhistsMasses)), std::ceil(sqrt(numhistsMasses)));
   for(int i=0; i<numhistsMasses; i++){
      canvasMass->cd(i+1);
      DataFitCompare(rootFileName, histnamesMasses.at(i), false);
   }

   TCanvas* canvasHeliAngles = new TCanvas("Angles helicity stystem","Angles helicity stystem",50, 70, 900,900);
   int numhistsHeliAngles = histnamesAnglesHeli.size();
   canvasHeliAngles->Divide(std::ceil(sqrt(numhistsHeliAngles)), std::ceil(sqrt(numhistsHeliAngles)));
   for(int i=0; i<numhistsHeliAngles; i++){
      canvasHeliAngles->cd(i+1);
      DataFitCompare(rootFileName, histnamesAnglesHeli.at(i), false);
   }

   TCanvas* canvasGJAngles = new TCanvas("Angles Gottfried Jackson stystem","Angles Gottfried Jackson stystem",80, 120,900,900);
   int numhistsGJAngles = histnamesAnglesGJ.size();
   canvasGJAngles->Divide(std::ceil(sqrt(numhistsGJAngles)), std::ceil(sqrt(numhistsGJAngles)));
   for(int i=0; i<numhistsGJAngles; i++){
      canvasGJAngles->cd(i+1);
      DataFitCompare(rootFileName, histnamesAnglesGJ.at(i), false);
   }

   if(saveImage){
      canvasMass->SaveAs("DataFitCompFileMasses.png");
      canvasHeliAngles->SaveAs("DataFitCompFileHeliAngles.png");
      canvasGJAngles->SaveAs("DataFitCompFileGJAngles.png");
   }
}

