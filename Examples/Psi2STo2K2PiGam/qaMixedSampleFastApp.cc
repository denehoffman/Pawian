#include <vector>
#include <map>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TMath.h"
#include <iostream>
#include "TCanvas.h"
#include "TF1.h"
#include <sstream>


TFile* fdata;
TNtuple* ntdata;
TFile* fpwamc;
TNtuple* ntpwamc;

TFile* fhistos;

TH1F* dataCosThetaKKpipi;
TH1F* pwamcCosThetaKKpipi;
TH1F* dataPhiKKpipi;
TH1F* pwamcPhiKKpipi;
TH1F* dataCosThetaKpipi;
TH1F* pwamcCosThetaKpipi;
TH1F* dataPhiKpipi;
TH1F* pwamcPhiKpipi;
TH1F* dataCosThetapipi;
TH1F* pwamcCosThetapipi;
TH1F* dataPhipipi;
TH1F* pwamcPhipipi;
TH1F* dataCosThetapi;
TH1F* pwamcCosThetapi;
TH1F* dataPhipi;
TH1F* pwamcPhipi;
TH1F* datainvMassKpipi;
TH1F* pwamcinvMassKpipi;
TH1F* datainvMasspipi;
TH1F* pwamcinvMasspipi;

TH1F* histofpull;

Int_t evtNo = 0;

float datacosThetaKKpipi;
float dataphiKKpipi;
float datacosThetaK1pi1pi2;
float dataphiK1pi1pi2;
float datacosThetapi1pi2;
float dataphipi1pi2;
float datacosThetapi1;
float dataphipi1;
float datainvmassK1pi1pi2;
float datainvmasspi1pi2;

float pwamccosThetaKKpipi;
float pwamcphiKKpipi;
float pwamccosThetaK1pi1pi2;
float pwamcphiK1pi1pi2;
float pwamccosThetapi1pi2;
float pwamcphipi1pi2;
float pwamccosThetapi1;
float pwamcphipi1;
float pwamcinvmassK1pi1pi2;
float pwamcinvmasspi1pi2;

float evtweight;


std::vector<float> veccosThetaKKpipi;
std::vector<float> vecphiKKpipi;
std::vector<float> veccosThetaK1pi1pi2;
std::vector<float> vecphiK1pi1pi2;
std::vector<float> veccosThetapi1pi2;
std::vector<float> vecphipi1pi2;
std::vector<float> veccosThetapi1;
std::vector<float> vecphipi1;
std::vector<float> vecinvmassK1pi1pi2;
std::vector<float> vecinvmasspi1pi2;
std::vector<float> vecevtweight;

std::vector<TH1F*> histVectData;
std::vector<TH1F*> histVectPwaMc;
std::vector<TH1F*> histVecthistofpull;

Int_t numOfEntriesData;
Int_t numOfEntriesPwaMc;

TString type = "mc";
bool drawFlag = true;
bool output = false;
bool longoutput = false;
bool silent = false;
TString method = "both";
using namespace std;

struct nextneighbours
{ 
  nextneighbours(  unsigned int theevtnumber1, unsigned int theevtnumber2, float thedistance, int thedatatype){
    evtnumber1 = theevtnumber1;
    evtnumber2 = theevtnumber2;
    distance = thedistance;
    datatype = thedatatype;
  }
 
  virtual ~nextneighbours(){};

  unsigned int evtnumber1;
  unsigned int evtnumber2;
  float distance;
  int datatype;
};

void fillHistograms();
void fillVectors(int start1, int start2, int end1, int end2);
void drawHistograms();
void chisq_method();
void mixed_sample();
float GetDistance(int entry1, int entry2);
bool Compare(const nextneighbours a, const nextneighbours b);


int main(){
  
  TString option="mixedsample,pwamc,nodraw,output";

  TString fileData="/data/sleipnir1/jansch/pwa_data/Psi2STo2K2PiGamPWA_data_2874events.root"; 
  //TString fileData="./Psi2STo2K2PiGamPWA_pwamc_FitParamHyp1_76695events.root";
  //TString fileData="./Psi2STo2K2PiGamPWA_pwamc_FitParamHyp9_43791events.root";
  //TString filePwaMc="./Psi2STo2K2PiGamPWA_pwamc_FitParamHyp1_76695events.root";
  //TString filePwaMc="./Psi2STo2K2PiGamPWA_pwamc_FitParamHyp9_43791events.root"; 
  //TString fileData="/data/sleipnir1/jansch/pwa_data/Psi2STo2K2PiGamPWA_pwamc_FitParamHyp1_739231events.root";
  TString filePwaMc="/data/sleipnir1/jansch/pwa_data/Psi2STo2K2PiGamPWA_pwamc_FitParamHyp1_739231events.root";
  cout << "--- ENTERING QAMIXEDSAMPLE" << endl;
  //   gROOT->SetStyle("Plain");
  //   gStyle->SetCanvasColor(0);
  //   gStyle->SetStatBorderSize(0);
  //   gStyle->SetPalette(1);
  //   gStyle->SetOptStat(1111);
  

  if (option.Contains("pwamc")) {type = "pwamc";}
  if (option.Contains("nodraw")) {drawFlag = false;}
  if (option.Contains("chisq")) {method = "chisq";}
  if (option.Contains("mixedsample")) {method = "mixedsample";}
  if (option.Contains("output")) {output = true;}
  if (option.Contains("longoutput")) {longoutput = true;}
  if (option.Contains("silent")) {silent = true;}

  fdata = new TFile(fileData,"READ");
  ntdata = (TNtuple*)fdata->Get("dataTuple");

  fpwamc = new TFile(filePwaMc,"READ");
  if(type=="pwamc") {
    cout << "Vergleich mit PWA-MC" << endl;
    ntpwamc = (TNtuple*)fpwamc->Get("dataTuple");
  }
  else if(type=="mc") {
    cout << "Vergleich mit MC" << endl;
    ntpwamc = (TNtuple*)fdata->Get("mcTuple");
  }
  cout << endl;

  numOfEntriesData = ntdata->GetEntries();
  cout << "Events in DATA Set: " << numOfEntriesData << endl;
  numOfEntriesPwaMc = ntpwamc->GetEntries();
  cout << "Events in PWAMC Set: " << numOfEntriesPwaMc << endl;

  dataCosThetaKKpipi = new TH1F("dataCosThetaKKpipi","dataCosThetaKKpipi",60,-1.,1.);
  dataPhiKKpipi = new TH1F("dataPhiKKpipi","dataPhiKKpipi",60,-TMath::Pi(),TMath::Pi());
  dataCosThetaKpipi = new TH1F("dataCosThetaKpipi","dataCosThetaKpipi",60,-1.,1.);
  dataPhiKpipi = new TH1F("dataPhiKpipi","dataPhiKpipi",60,-TMath::Pi(),TMath::Pi());
  dataCosThetapipi = new TH1F("dataCosThetapipi","dataCosThetapipi",60,-1.,1.);
  dataPhipipi = new TH1F("dataPhipipi","dataPhipipi",60,-TMath::Pi(),TMath::Pi());
  dataCosThetapi = new TH1F("dataCosThetapi","dataCosThetapi",60,0.,1.);
  dataPhipi = new TH1F("dataPhipi","dataPhipi",60,0.,TMath::Pi());
  datainvMassKpipi = new TH1F("datainvMassKpipi", "datainvMassKpipi", 60, 0.7, 3.0);
  datainvMasspipi = new TH1F("datainvMasspipi", "datainvMasspipi", 60, 0.2, 2.8);
   
  pwamcCosThetaKKpipi = new TH1F("pwamcCosThetaKKpipi","pwamcCosThetaKKpipi",60,-1.,1.);
  pwamcPhiKKpipi = new TH1F("pwamcPhiKKpipi","pwamcPhiKKpipi",60,-TMath::Pi(),TMath::Pi());
  pwamcCosThetaKpipi = new TH1F("pwamcCosThetaKpipi","pwamcCosThetaKpipi",60,-1.,1.);
  pwamcPhiKpipi = new TH1F("pwamcPhiKpipi","pwamcPhiKpipi",60,-TMath::Pi(),TMath::Pi());
  pwamcCosThetapipi = new TH1F("pwamcCosThetapipi","pwamcCosThetapipi",60,-1.,1.);
  pwamcPhipipi = new TH1F("pwamcPhipipi","pwamcPhipipi",60,-TMath::Pi(),TMath::Pi());
  pwamcCosThetapi = new TH1F("pwamcCosThetapi","pwamcCosThetapi",60,0.,1.);
  pwamcPhipi = new TH1F("pwamcPhipi","pwamcPhipi",60,0.,TMath::Pi());
  pwamcinvMassKpipi = new TH1F("pwamcinvMassKpipi", "pwamcinvMassKpipi", 60, 0.7, 3.0);
  pwamcinvMasspipi = new TH1F("pwamcinvMasspipi", "pwamcinvMasspipi", 60, 0.2, 2.8);
  
  histVectData.push_back(dataCosThetaKKpipi);
  histVectData.push_back(dataPhiKKpipi);
  histVectData.push_back(dataCosThetaKpipi);
  histVectData.push_back(dataPhiKpipi);
  histVectData.push_back(dataCosThetapipi);
  histVectData.push_back(dataPhipipi);
  histVectData.push_back(dataCosThetapi);
  histVectData.push_back(dataPhipi);
  histVectData.push_back(datainvMassKpipi);
  histVectData.push_back(datainvMasspipi);

  histVectPwaMc.push_back(pwamcCosThetaKKpipi);
  histVectPwaMc.push_back(pwamcPhiKKpipi);
  histVectPwaMc.push_back(pwamcCosThetaKpipi);
  histVectPwaMc.push_back(pwamcPhiKpipi);
  histVectPwaMc.push_back(pwamcCosThetapipi);
  histVectPwaMc.push_back(pwamcPhipipi);
  histVectPwaMc.push_back(pwamcCosThetapi);
  histVectPwaMc.push_back(pwamcPhipi);
  histVectPwaMc.push_back(pwamcinvMassKpipi);
  histVectPwaMc.push_back(pwamcinvMasspipi);

  fillHistograms();

  if(drawFlag) {
    cout << "--- Histogramme gezeichnet." << endl;
    drawHistograms();
    cout << endl;
    cout << endl;
    cout << "-----------------------------------------" << endl;
  }
  if(method=="chisq" || method=="both") {
    cout << "--- CHISQ-METHOD-ERGEBNISSE:" << endl;
    cout << endl;
    chisq_method();
    cout << "-----------------------------------------" << endl;
    cout << endl;
  }
  if(method=="mixedsample" || method=="both") {
    cout << "--- MIXED-SAMPLE-ERGEBNISSE:" << endl;
    cout << endl;
    mixed_sample();
    cout << "-----------------------------------------" << endl;
    cout << endl;
  }
  cout << "--- EXITING QAMIXEDSAMPLE" << endl;

  return 1;

}


bool Compare(const nextneighbours a, const nextneighbours b) {
  return a.distance < b.distance;
}


void fillVectors(int start1, int end1, int start2, int end2) {

  veccosThetaKKpipi.clear();
  vecphiKKpipi.clear();
  veccosThetaK1pi1pi2.clear();
  vecphiK1pi1pi2.clear();
  veccosThetapi1pi2.clear();
  vecphipi1pi2.clear();
  veccosThetapi1.clear();
  vecphipi1.clear();
  vecinvmassK1pi1pi2.clear();
  vecinvmasspi1pi2.clear();
  cout << "Vectors cleared." << endl;
  cout << endl;
  cout << "Read branches vom ntdata." << endl;

  ntdata->SetBranchAddress("costhetaKKpipi", &datacosThetaKKpipi);
  ntdata->SetBranchAddress("phiKKpipi", &dataphiKKpipi);
  ntdata->SetBranchAddress("costhetaK1pi1pi2", &datacosThetaK1pi1pi2);
  ntdata->SetBranchAddress("phiK1pi1pi2", &dataphiK1pi1pi2);
  ntdata->SetBranchAddress("costhetapipiViaK1pipi", &datacosThetapi1pi2);
  ntdata->SetBranchAddress("phipipiViaK1pipi", &dataphipi1pi2);
  ntdata->SetBranchAddress("costhetapiViapipi", &datacosThetapi1);
  ntdata->SetBranchAddress("phipiViapipi", &dataphipi1);
  ntdata->SetBranchAddress("mk1pi1pi2", &datainvmassK1pi1pi2);
  ntdata->SetBranchAddress("mpipi", &datainvmasspi1pi2);
  ntpwamc->SetBranchAddress("costhetaKKpipi", &pwamccosThetaKKpipi);
  ntpwamc->SetBranchAddress("phiKKpipi", &pwamcphiKKpipi);
  ntpwamc->SetBranchAddress("costhetaK1pi1pi2", &pwamccosThetaK1pi1pi2);
  ntpwamc->SetBranchAddress("phiK1pi1pi2", &pwamcphiK1pi1pi2);
  ntpwamc->SetBranchAddress("costhetapipiViaK1pipi", &pwamccosThetapi1pi2);
  ntpwamc->SetBranchAddress("phipipiViaK1pipi", &pwamcphipi1pi2);
  ntpwamc->SetBranchAddress("costhetapiViapipi", &pwamccosThetapi1);
  ntpwamc->SetBranchAddress("phipiViapipi", &pwamcphipi1);
  ntpwamc->SetBranchAddress("mk1pi1pi2", &pwamcinvmassK1pi1pi2);
  ntpwamc->SetBranchAddress("mpipi", &pwamcinvmasspi1pi2);
  ntpwamc->SetBranchAddress("weight",&evtweight);

  for(int i = start1; i<end1; i++) {
    ntdata->GetEntry(i);
    veccosThetaKKpipi.push_back(datacosThetaKKpipi);
    vecphiKKpipi.push_back(dataphiKKpipi);
    veccosThetaK1pi1pi2.push_back(datacosThetaK1pi1pi2);
    vecphiK1pi1pi2.push_back(dataphiK1pi1pi2);
    veccosThetapi1pi2.push_back(datacosThetapi1pi2);
    vecphipi1pi2.push_back(dataphipi1pi2);
    veccosThetapi1.push_back(datacosThetapi1);
    vecphipi1.push_back(dataphipi1);
    vecinvmassK1pi1pi2.push_back(datainvmassK1pi1pi2);
    vecinvmasspi1pi2.push_back(datainvmasspi1pi2);
  }
  cout << "Data Events filled to vectors." << endl;

  cout << "Read branches vom ntpwamc." << endl;

  for(int i = start2; i<end2; i++) {
    ntpwamc->GetEntry(i);
    veccosThetaKKpipi.push_back(pwamccosThetaKKpipi);
    vecphiKKpipi.push_back(pwamcphiKKpipi);
    veccosThetaK1pi1pi2.push_back(pwamccosThetaK1pi1pi2);
    vecphiK1pi1pi2.push_back(pwamcphiK1pi1pi2);
    veccosThetapi1pi2.push_back(pwamccosThetapi1pi2);
    vecphipi1pi2.push_back(pwamcphipi1pi2);
    veccosThetapi1.push_back(pwamccosThetapi1);
    vecphipi1.push_back(pwamcphipi1);
    vecinvmassK1pi1pi2.push_back(pwamcinvmassK1pi1pi2);
    vecinvmasspi1pi2.push_back(pwamcinvmasspi1pi2);
  }
  cout << "PwaMC Events filled to vectors." << endl;
  cout << "Filled " << veccosThetaKKpipi.size() << " events in total." << endl;

}



void mixed_sample() {

  if(!silent) cout << "Entering mixed-sample method" << endl;
  if(!silent) cout << endl;

  std::vector<int> vec_limit_data;
  //vec_limit_data.push_back(100);
  //vec_limit_data.push_back(200);
  //vec_limit_data.push_back(500);
  vec_limit_data.push_back(1000);
  //vec_limit_data.push_back(2000);
  //vec_limit_data.push_back(2874);
  //vec_limit_data.push_back(3000);
  //vec_limit_data.push_back(5000);
  //vec_limit_data.push_back(7500);
  // vec_limit_data.push_back(10000);
  
  std::vector<int> vec_limit_pwamc;
  //vec_limit_pwamc.push_back(1000);
  //vec_limit_pwamc.push_back(2000);
  //vec_limit_pwamc.push_back(3000);
  //vec_limit_pwamc.push_back(4000);
  //vec_limit_pwamc.push_back(5000);
  vec_limit_pwamc.push_back(7500);
  vec_limit_pwamc.push_back(10000);
  vec_limit_pwamc.push_back(15000);
  vec_limit_pwamc.push_back(20000);

  for(unsigned int i = 0; i<vec_limit_data.size(); i++) {
    cout << "Cycle " << i << " will have " << vec_limit_data[i] << " data events." << endl;
  }

  for(unsigned int k = 0; k<vec_limit_pwamc.size(); k++) {

    //unsigned int limit_data = 2874;   // numOfEntriesData (usually 2874);
    unsigned int limit_data = 1000; // vec_limit_data[k];
    unsigned int limit_pwamc = vec_limit_pwamc[k]; // 10*limit_data;  // numOfEntriesPwaMc;
    string sourcename1 = "data";
    string sourcename2 = "Hyp1pwamc";
    unsigned int limit_runs_lower = 1;
    unsigned int limit_runs_upper = numOfEntriesPwaMc/limit_pwamc;
    unsigned int numberofneighbours = 10;
    std::vector<float> fpullvector;

    cout << "Cycle: " << k << endl;
    cout << "Using " << limit_data << " DATA events per run." << endl;
    cout << "Using " << limit_pwamc << " PWAMC events per run." << endl;
    cout << "Performing " << limit_runs_upper - limit_runs_lower << " runs." << endl;
       
    
    for(unsigned int run = limit_runs_lower; run < limit_runs_upper; run++) {

      cout << "Starting run " << run << endl;
      if(!silent) cout << endl;
      
      float distance = 1E6;
      int type = 0;
      float T=0;
      nextneighbours maxNeighbours(0,0,1000000.,0);
      nextneighbours tmpNeighbours(0,0,100000,99);
      std::vector<struct nextneighbours> neighbourList;
      unsigned int numOfEntriesTotalUsed = limit_data + limit_pwamc;

      fillVectors(0, limit_data, run*limit_pwamc, (run+1)*limit_pwamc);
      cout << "Vectors have been filled." << endl;
      
      if(!silent) {
	cout << endl;
	cout << "Number of Entries in DATA: " << numOfEntriesData << " events available (used " << limit_data << ")"<< endl;
	cout << "Number of Entries in PWAMC: " << numOfEntriesPwaMc << " events available (used " << limit_pwamc << " starting at " << run*limit_pwamc << " to " << ((run+1)*limit_pwamc)-1 << ")" << endl;
	cout << "Total Number of used Entries: " << numOfEntriesTotalUsed << endl;
	cout << endl;
      }

      for(unsigned int i = 0; i < numOfEntriesTotalUsed; i++) {
        neighbourList.clear();
	for(unsigned int nn = 0; nn < numberofneighbours; nn++) {
	  neighbourList.push_back(maxNeighbours);
	}
	for(unsigned int j = 0; j < numOfEntriesTotalUsed; j++) {
          if(i!=j) {
	    std::sort(neighbourList.begin(),neighbourList.end(), Compare);
	    if(i<limit_data && j<limit_data)
	      {
		distance = GetDistance(i,j);
		type = 1;
		tmpNeighbours.evtnumber1 = i;
		tmpNeighbours.evtnumber2 = j;
		tmpNeighbours.distance = distance;
		tmpNeighbours.datatype = type;
		}
	    else if (i<limit_data && j>=limit_data) 
	      {
		distance = GetDistance(i,j);
                type = 0;
		tmpNeighbours.evtnumber1 = i;
                tmpNeighbours.evtnumber2 = j;
                tmpNeighbours.distance = distance;
                tmpNeighbours.datatype = type;
	      }
	    else if (i>=limit_data && j<limit_data) 
	      {
		distance = GetDistance(i,j);
                type = 0;
		tmpNeighbours.evtnumber1 = i;
                tmpNeighbours.evtnumber2 = j;
                tmpNeighbours.distance = distance;
                tmpNeighbours.datatype = type;
	      }
	    else if (i>=limit_data && j>=limit_data)
              {
                distance = GetDistance(i,j);
                type = 1;
		tmpNeighbours.evtnumber1 = i;
                tmpNeighbours.evtnumber2 = j;
                tmpNeighbours.distance = distance;
                tmpNeighbours.datatype = type;
	      }
            if(distance < neighbourList[numberofneighbours-1].distance) { 
	      neighbourList[numberofneighbours-1] = tmpNeighbours;
	    }
          }
        }
	std::sort(neighbourList.begin(),neighbourList.end(), Compare);
	if(longoutput) cout << neighbourList.size() << " next neighbours of event " << i << " are" << endl;
	for(unsigned int k = 0; k < neighbourList.size(); k++) {
	  if(longoutput) cout << neighbourList[k].evtnumber1 << "\t" << neighbourList[k].evtnumber2 << "\t" << neighbourList[k].distance << " \t\t" << neighbourList[k].datatype << endl;
	  T=T+neighbourList[k].datatype;
	}
	neighbourList.clear();
	if(((i+1) % 1000) == 0 && !silent) cout << "Finished " << i+1 << "th event." << endl;
      }

      float Tnorm = T/(float(numberofneighbours)*float(numOfEntriesTotalUsed));
      float mu = ((float(limit_data)*(float(limit_data)-1)) + (float(limit_pwamc)*(float(limit_pwamc)-1))) / (float(numOfEntriesTotalUsed)*(float(numOfEntriesTotalUsed)-1));
      float variancesquared = 
	1 / (float(numOfEntriesTotalUsed)*float(numberofneighbours)) 
	* 
	( 
	 (float(limit_data)*float(limit_pwamc))/(pow(float(numOfEntriesTotalUsed), 2)) 
	 + 
	 4*(pow(float(limit_data), 2)*pow(float(limit_pwamc), 2))/(pow(float(numOfEntriesTotalUsed),4)) 
	 );
      float fpull = (Tnorm-mu)/sqrt(variancesquared);

      if(!silent) cout << "T: " << T << endl;
      if(!silent) cout << "Tnorm: " << Tnorm << endl;
      if(!silent) cout << "Erwartungswert mu: " << mu << endl;
      if(!silent) {
	cout << "Varianz: " << variancesquared << endl;
	cout << endl; 
	cout << "Differenz zwischen Resultat und Erwartungswert: " << Tnorm-mu << endl;
      }
      cout << "f_pull in run " << run << ": " << fpull << endl;
      cout << endl;
      fpullvector.push_back(fpull);
    }

    std::stringstream out;
    out << limit_data;
    string limit_data_string = out.str();
    out.str("");
    out.clear();
    out << limit_pwamc;
    string limit_pwamc_string = out.str();
    out.str("");
    out.clear();
    out << limit_runs_upper - limit_runs_lower;
    string runs_string = out.str();
    
    TString histname = "histofpull_"+limit_data_string+sourcename1+"_"+limit_pwamc_string+sourcename2+"_"+runs_string+"runs.pdf";
    //cout << histname << endl;
    histofpull = new TH1F("histofpull", histname, 100, -10, 10);

    cout << "--------------------------" << endl;
    cout << fpullvector.size() << " fpulls have been calculated." << endl;
    for(unsigned int i = 0; i < fpullvector.size(); i++) {
      cout << "Calculated fpull:\t" << fpullvector[i] << endl;
      histofpull->Fill( fpullvector[i] );
    }
        
    TCanvas* chist1 = new TCanvas("chist1","chist1",1200,800);
    histofpull->SetLineColor(2);
    histofpull->Draw();
    
    TF1* fit1 = new TF1("fit1","gaus(0)",-10.,10.);
    fit1->SetParameters(5,histofpull->GetMean(),histofpull->GetRMS());
    fit1->SetLineColor(4);
    fit1->SetLineWidth(3);
    histofpull->Fit("fit1","LMRS","",-10.,10.);

    chist1->Print(histname);

    histVecthistofpull.push_back(histofpull);

    TString rootfile = "histofpull_"+limit_data_string+sourcename1+"_"+limit_pwamc_string+sourcename2+"_"+runs_string+"runs.root";
    fhistos = new TFile(rootfile,"RECREATE");
    for(unsigned int i = 0; i < histVecthistofpull.size(); i++) 
      {
	histVecthistofpull[i]->Write();
      }

    histofpull->Reset();
    fpullvector.clear();
    
        
  }

}




float GetDistance(int entry1, int entry2) {
  
  float distance = sqrt( 
			 pow((veccosThetaKKpipi[entry1]-veccosThetaKKpipi[entry2])/2.,2) +
			 pow((vecphiKKpipi[entry1]-vecphiKKpipi[entry2])/(2*3.1415),2) +
			 pow((veccosThetaK1pi1pi2[entry1]-veccosThetaK1pi1pi2[entry2])/2. ,2) +
			 pow((vecphiK1pi1pi2[entry1]-vecphiK1pi1pi2[entry2])/(2*3.1415) ,2) +
			 pow((veccosThetapi1pi2[entry1]-veccosThetapi1pi2[entry2])/2. ,2) +
			 pow((vecphipi1pi2[entry1]-vecphipi1pi2[entry2])/(2*3.1415) ,2) +
			 pow((veccosThetapi1[entry1]-veccosThetapi1[entry2]) ,2) +
			 pow((vecphipi1[entry1]-vecphipi1[entry2])/(3.1415) ,2) +
			 pow((vecinvmassK1pi1pi2[entry1]-vecinvmassK1pi1pi2[entry2] )/2.2 ,2) +
			 pow((vecinvmasspi1pi2[entry1]-vecinvmasspi1pi2[entry2])/2.3 ,2)
			 );
  
  return distance;

}



void chisq_method() {

  cout << "Entering chisq-method" << endl;

  Double_t chisq = 0;

  Int_t constraints = 50; // 12-01-06_FitParamHyp8_base_K0K2_Hyp5/rerun_correctedampK0K2/Psi2STo2K2PiGamPWA_noK1_1680Hyp7.root
  // float constraints = 94;   // 12-01-06_FitParamHyp8_base_K0K2_Hyp5/rerun_correctedampK0K2/Psi2STo2K2PiGamPWA_noK1_1680Hyp7_noK1_2300Hyp7.root
  Int_t ndf;
  int numberofbins = 0;
  
  for(unsigned int i=0; i<histVectData.size(); i++) {
    int currentNoOfBins=0;

    for (int j=0; j<(histVectData[i]->GetSize()); j++) {
      float binContentData=float(histVectData[i]->GetBinContent(j));
      if (binContentData>0.){
	chisq += pow( binContentData - (histVectPwaMc[i]->GetBinContent(j)),2)/binContentData;
	currentNoOfBins++;
      }
      // cout << "chi^2 after bin " << j << " of histogram " << histVectData[i]->GetName() << ": " << chisq << endl;
    }
    cout << "chi^2 after " << histVectData[i]->GetName() << ": " << chisq << endl;
    
    numberofbins += currentNoOfBins;
  }   

    ndf = numberofbins - constraints;

    cout << "Total chi^2:\t" << chisq << endl;    
    cout << "Total number of bins:\t" << numberofbins << endl;
    cout << "Number of constraints:\t" << constraints << endl;
    cout << "Degrees of freedom:\t" << ndf << endl;
    cout << "chi^2/ndf:\t" << chisq/double(ndf) << endl;
    cout << "CL:\t" << TMath::Prob(chisq, ndf) << endl;
}






void fillHistograms() {

  Int_t maxEvtNo = 0;
  for(Int_t i=0; i<max(numOfEntriesPwaMc,numOfEntriesData); i++)
    {
      if(i<numOfEntriesData){
      	ntdata->GetEntry(i);

	dataCosThetaKKpipi->Fill(datacosThetaKKpipi);
        dataPhiKKpipi->Fill(dataphiKKpipi);	

        dataCosThetaKpipi->Fill(datacosThetaK1pi1pi2);	
        dataPhiKpipi->Fill(dataphiK1pi1pi2);	

	dataCosThetapipi->Fill(datacosThetapi1pi2);
       	dataPhipipi->Fill(dataphipi1pi2);
 
        dataCosThetapi->Fill(fabs(datacosThetapi1));                                              	
        dataPhipi->Fill(fabs(dataphipi1));                                              	
   
        datainvMassKpipi->Fill(datainvmassK1pi1pi2);
        datainvMasspipi->Fill(datainvmasspi1pi2);
      }

      ntpwamc->GetEntry(i);

      pwamcCosThetaKKpipi->Fill(pwamccosThetaKKpipi, evtweight);
      pwamcPhiKKpipi->Fill(pwamcphiKKpipi, evtweight);	
      
      pwamcCosThetaKpipi->Fill(pwamccosThetaK1pi1pi2, evtweight);	
      pwamcPhiKpipi->Fill(pwamcphiK1pi1pi2, evtweight);	
      
      pwamcCosThetapipi->Fill(pwamccosThetapi1pi2, evtweight);
      pwamcPhipipi->Fill(pwamcphipi1pi2, evtweight);
      
      pwamcCosThetapi->Fill(fabs(pwamccosThetapi1), evtweight);                                              	
      pwamcPhipi->Fill(fabs(pwamcphipi1), evtweight);                                              	
      
      pwamcinvMassKpipi->Fill(pwamcinvmassK1pi1pi2, evtweight);
      pwamcinvMasspipi->Fill(pwamcinvmasspi1pi2, evtweight);
      
      if(maxEvtNo<evtNo) maxEvtNo=evtNo;
    }


}




void drawHistograms() {

  double scalefactor = double(ntdata->GetEntries())/double(ntpwamc->GetEntries()); 

  Int_t rebin = 1;

  TCanvas* cmain = new TCanvas("cmain","cmain",1400,600);
  cmain->Divide(5,2);
  for(unsigned int i=0; i<histVectData.size(); i++) {
    cmain->cd(i+1);
    histVectData[i]->Rebin(rebin);
    histVectData[i]->SetLineColor(2);
    histVectData[i]->Draw("E");
    histVectData[i]->SetMinimum(0);
    histVectPwaMc[i]->Rebin(rebin);
    histVectPwaMc[i]->Scale(scalefactor);
    histVectPwaMc[i]->Draw("same");
  } 

}



