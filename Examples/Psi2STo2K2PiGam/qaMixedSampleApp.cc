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

TFile* fdata;
TNtuple* ntdata;
TFile* fpwamc;
TNtuple* ntpwamc;

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

std::vector<TH1F*> histVectData;
std::vector<TH1F*> histVectPwaMc;

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
  nextneighbours(  unsigned int theevtnumber1, unsigned int theevtnumber2, float thedistance, bool thedatatype){
    evtnumber1 = theevtnumber1;
    evtnumber2 = theevtnumber2;
    distance = thedistance;
    datatype = thedatatype;
  }
 
  virtual ~nextneighbours(){};

  unsigned int evtnumber1;
  unsigned int evtnumber2;
  float distance;
  bool datatype;
};

void fillHistograms();
void drawHistograms();
void chisq_method();
void mixed_sample();
float GetDistanceDataToData(int dataentry1, int dataentry2);
float GetDistanceDataToPwaMc(int dataentry, int pwamcentry);
float GetDistancePwaMcToPwaMc(int pwamcentry1, int pwamcentry2);
bool Compare(const nextneighbours a, const nextneighbours b);

/*
void initNtData(){

  if (0 != fdata){
     delete fdata;
    fdata=0;
  }

  TString fileData="./Psi2STo2K2PiGamPWA_data.root"; 
  fdata = new TFile(fileData,"READ");
  ntdata = (TNtuple*)fdata->Get("dataTuple");
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
}


void initNtPwamc(){

  if (0 != fpwamc){
     delete fpwamc;
    fpwamc=0;
  }

  TString filePwaMc="./Psi2STo2K2PiGamPWA_pwamc__FitParamHyp8_base_K0K2_newflatte__10000events.root";
  fpwamc = new TFile(filePwaMc,"READ");
  ntpwamc = (TNtuple*)fpwamc->Get("dataTuple");
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
}
*/

int main(){
  
  TString fileData="./Psi2STo2K2PiGamPWA_data.root"; 
  //TString fileData="./Psi2STo2K2PiGamPWA_pwamc__FitParamHyp8_base_K0K2_newflatte__10000events.root";
  TString filePwaMc="./Psi2STo2K2PiGamPWA_pwamc__FitParamHyp8_base_K0K2_newflatte__10000events.root"; 
  TString option="mixedsample,pwamc,nodraw,silent";

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

  numOfEntriesData = ntdata->GetEntries();
  numOfEntriesPwaMc = ntpwamc->GetEntries();


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
  
  histofpull = new TH1F("histofpull", "histofpull", 100, -10, 10);

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


void mixed_sample() {

  if(!silent) cout << "Entering mixed-sample method" << endl;
  if(!silent) cout << endl;

  unsigned int limit_runs_lower = 0;
  unsigned int limit_runs_upper = 10;
  cout << "Performing " << limit_runs_upper - limit_runs_lower << " runs." << endl;
  std::vector<float> fpullvector;
  unsigned int limit_data = 1000;   // numOfEntriesData;
  unsigned int limit_pwamc = 1000;  // numOfEntriesPwaMc;
  unsigned int numberofneighbours = 10;


  for(unsigned int run = limit_runs_lower; run < limit_runs_upper; ++run) {

    std::vector<struct nextneighbours> neighbourList;
    cout << "Starting run " << run+1 << endl;

    if(!silent) cout << endl;

    float distance = 1E6;
    float T=0;
    // unsigned int limit_data = 2500;   // numOfEntriesData;
    // unsigned int limit_pwamc = 5000;  // numOfEntriesPwaMc;
    // unsigned int numberofneighbours = 10;
    unsigned int pwaoffset = 100000;

    if(!silent) cout << "--- Calculating distances DATA <-> DATA" << endl;
    for(unsigned int i = 0; i<limit_data; i++) {
      for(unsigned int j = i; j<limit_data; j++) {
        //cout << run << ", loop 1:\t" << i << "\t" << j << endl;
	if(i!=j) {
	  distance = GetDistanceDataToData(i,j);
          nextneighbours tmpNeighbours(i,j,distance,true);
	  neighbourList.push_back(tmpNeighbours);
	  if(neighbourList.size() % 500000 == 0) cout << neighbourList.size() << "th distance has been calculated." << endl;
	}
      }
    }
    
    if(!silent)  cout << "--- Calculating distances DATA <-> PWAMC" << endl;
    for(unsigned int i = 0; i<limit_data; i++) {
      for(unsigned int j = 0; j<limit_pwamc; j++) {
        //cout <<  run << ", loop 2:\t" << i << "\t" << j << endl;
  	distance = GetDistanceDataToPwaMc(i,j+run*limit_pwamc);
	nextneighbours tmpNeighbours(i,j+pwaoffset,distance,false);
	neighbourList.push_back(tmpNeighbours);
	// neighbourList.push_back((struct nextneighbours) {i,j+pwaoffset,distance,false});
        if(neighbourList.size() % 500000 == 0) cout << neighbourList.size() << "th distance has been calculated." << endl;
      }
    }
    
    if(!silent) cout << "--- Calculating distances PWAMC <-> PWAMC" << endl;
    for(unsigned int i = 0; i<limit_pwamc; i++) {
      for(unsigned int j = i; j<limit_pwamc; j++) {
        //cout << run << ", loop 3:\t" << i << "\t" << j << endl;
	if(i!=j) {
	  distance = GetDistancePwaMcToPwaMc(i+run*limit_pwamc,j+run*limit_pwamc);
	  nextneighbours tmpNeighbours(i+pwaoffset,j+pwaoffset,distance,true);
	  neighbourList.push_back(tmpNeighbours);
	  // neighbourList.push_back((struct nextneighbours) {i+pwaoffset,j+pwaoffset,distance,true});
	  if(neighbourList.size() % 500000 == 0) cout << neighbourList.size() << "th distance has been calculated." << endl;
	}
      }
    }
    
    if(!silent) {
      cout << endl;
      cout << "How many distances have been calculated?\t" << neighbourList.size() << endl;
      cout << endl;
      cout << endl;
    }
    
    if(longoutput) {
      cout << "--- Distances:" << endl;
      for(unsigned int i = 0; i < neighbourList.size(); i++) {
	cout << neighbourList[i].evtnumber1 << "\t" << neighbourList[i].evtnumber2 << "\t" << neighbourList[i].distance << " \t\t" << neighbourList[i].datatype << endl;
      }
    }
    if(!silent) cout << endl;
    
    std::sort(neighbourList.begin(),neighbourList.end(), Compare);
    if(!silent) cout << "neighbourList has been sorted!" << endl;
    if(!silent) cout << endl;
    
    if(longoutput) {
      cout << "--- Distances after sorting:" << endl;
      for(unsigned int i = 0; i < neighbourList.size(); i++) {
	cout << neighbourList[i].evtnumber1 << "\t" << neighbourList[i].evtnumber2 << "\t" << neighbourList[i].distance << " \t\t" << neighbourList[i].datatype << endl;
      }
    }
    if(!silent)  cout << endl;
    
    float numOfEntriesTotal = float(limit_data + limit_pwamc);
    
    unsigned int m;
    
    for(unsigned int i = 0; i < limit_data; i++) {
      if(output) cout << "--- Checking neighbours of event " << i << endl;
      m = 0;
      for(unsigned int j = 0; j < neighbourList.size(); j++) 
	{
	  if (m >= numberofneighbours) continue;
	  if(longoutput) cout << "Checking Event in neighbourList for data: " << j << endl;
	  if(neighbourList[j].evtnumber1 == i || neighbourList[j].evtnumber2 == i) 
	    {
	      m++;
	      if(longoutput) cout << j << "th event is the " << m << "th neighbour of event " << i;
	      if(neighbourList[j].datatype==true) 
		{
		  if(longoutput) cout << " and belongs to the same sample! Adding +1 to T" << endl;
		  T=T+1.;
		}
	      else 
		{
		  if(longoutput) cout << ", but doesn't belong to the same sample." << endl;
		}
	    }
	  else 
	    {
	      if(longoutput) cout << "Not a neighbour of event " << i << endl;
	    }
	}
      if(output) cout << "Event " << i << " has been analysed with " << m << " neighbours."  << endl;
      if(output) cout << "T: " << T << endl;
      if(output) cout << endl;
    }

    for(unsigned int i = 0; i < limit_pwamc; i++) {
      if(output) cout << "--- Checking neighbours of event " << i+pwaoffset << endl;
      m = 0;
      
      for(unsigned int j = 0; j < neighbourList.size(); j++) 
	{
	  if (m >= numberofneighbours) continue;
	  if(longoutput) cout << "Checking Event in neighbourList for pwamc: " << j << endl;
	  if(neighbourList[j].evtnumber1 == i+pwaoffset || neighbourList[j].evtnumber2 == i+pwaoffset) 
	    {
	      m++;
	      if(longoutput) cout << j << "th event is the " << m << "th neighbour of event " << i+pwaoffset;
	      if(neighbourList[j].datatype==true) 
		{
		  if(longoutput) cout << " and belongs to the same sample! Adding +1 to T" << endl;
		  T=T+1;
		}
	      else 
		{
		  if(longoutput) cout << ", but doesn't belong to the same sample." << endl;
		}
	    }
	  else 
	    {
	      if(longoutput) cout << "Not a neighbour of event " << i+pwaoffset << endl;
	    }

	}

      if(output) cout << "Event " << i+pwaoffset << " has been analysed." << endl;
      if(output) cout << "T: " << T << endl;
      if(output) cout << endl;
    }


    if(!silent) {
      cout << endl;
      cout << "Number of Entries in DATA: " << numOfEntriesData << " events available (used " << limit_data << ")"<< endl;
      cout << "Number of Entries in PWAMC: " << numOfEntriesPwaMc << " events available (used " << limit_pwamc << ")"<< endl;
      cout << "Total Number of Entries: " << numOfEntriesTotal << endl;
      cout << endl;
    }
    
    numOfEntriesData = limit_data;
    numOfEntriesPwaMc = limit_pwamc;
    
    float Tnorm = 0;
    Tnorm = T/(float(numberofneighbours)*float(numOfEntriesTotal));
    if(!silent) cout << "Tnorm: " << Tnorm << endl;
    
    
    float mu = ((float(numOfEntriesData)*(float(numOfEntriesData)-1)) + (float(numOfEntriesPwaMc)*(float(numOfEntriesPwaMc)-1))) / (float(numOfEntriesTotal)*(float(numOfEntriesTotal)-1));
    if(!silent) cout << "Erwartungswert mu: " << mu << endl;
    
    
    float variancesquared = 
      1 / (float(numOfEntriesTotal)*float(numberofneighbours)) 
      * 
      ( 
       (float(numOfEntriesData)*float(numOfEntriesPwaMc))/(pow(float(numOfEntriesTotal), 2)) 
       + 
       4*(pow(float(numOfEntriesData), 2)*pow(float(numOfEntriesPwaMc), 2))/(pow(float(numOfEntriesTotal),4)) 
       );

    float fpull = (Tnorm-mu)/sqrt(variancesquared);
    
    if(!silent) {
      cout << "Varianz: " << variancesquared << endl;
      cout << endl; 
      cout << "Differenz zwischen Resultat und Erwartungswert: " << Tnorm-mu << endl;
    }
    cout << "f_pull in run " << run+1 << ": " << fpull << endl;
    cout << endl;
    fpullvector.push_back(fpull);
  }
  
  cout << "--------------------------" << endl;
  cout << fpullvector.size() << endl;
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
  chist1->Print("histofpull.pdf");

}




float GetDistanceDataToData(int dataentry1, int dataentry2) {
  //  cout << "Distance from Data Event " << dataentry1 << " to Data Event " << dataentry2 << endl;
  
  ntdata->GetEntry(dataentry1);
 float datacosThetaKKpipi_1 = datacosThetaKKpipi;
  float dataphiKKpipi_1 = dataphiKKpipi;
  float datacosThetaK1pi1pi2_1 = datacosThetaK1pi1pi2;
  float dataphiK1pi1pi2_1 = dataphiK1pi1pi2;
  float datacosThetapi1pi2_1 = datacosThetapi1pi2;
  float dataphipi1pi2_1 = dataphipi1pi2;
  float datacosThetapi1_1 = datacosThetapi1;
  float dataphipi1_1 = dataphipi1;
  float datainvmassK1pi1pi2_1 = datainvmassK1pi1pi2;
  float datainvmasspi1pi2_1 = datainvmasspi1pi2;
  ntdata->GetEntry(dataentry2);
  float datacosThetaKKpipi_2 = datacosThetaKKpipi;
  float dataphiKKpipi_2 = dataphiKKpipi;
  float datacosThetaK1pi1pi2_2 = datacosThetaK1pi1pi2;
  float dataphiK1pi1pi2_2 = dataphiK1pi1pi2;
  float datacosThetapi1pi2_2 = datacosThetapi1pi2;
  float dataphipi1pi2_2 = dataphipi1pi2;
  float datacosThetapi1_2 = datacosThetapi1;
  float dataphipi1_2 = dataphipi1;
  float datainvmassK1pi1pi2_2 = datainvmassK1pi1pi2;
  float datainvmasspi1pi2_2 = datainvmasspi1pi2;
  
  float distance = sqrt( 
			 pow((datacosThetaKKpipi_1-datacosThetaKKpipi_2)/2.,2) +
			 pow((dataphiKKpipi_1-dataphiKKpipi_2)/(2*3.1415),2) +
			 pow((datacosThetaK1pi1pi2_1-datacosThetaK1pi1pi2_2)/2. ,2) +
			 pow((dataphiK1pi1pi2_1-dataphiK1pi1pi2_2)/(2*3.1415) ,2) +
			 pow((datacosThetapi1pi2_1-datacosThetapi1pi2_2)/2. ,2) +
			 pow((dataphipi1pi2_1-dataphipi1pi2_2)/(2*3.1415) ,2) +
			 pow((datacosThetapi1_1-datacosThetapi1_2) ,2) +
			 pow((dataphipi1_1-dataphipi1_2)/(3.1415) ,2) +
			 pow((datainvmassK1pi1pi2_1-datainvmassK1pi1pi2_2 )/2.2 ,2) +
			 pow((datainvmasspi1pi2_1-datainvmasspi1pi2_2)/2.3 ,2)
			 );
  
  return distance;

}



float GetDistanceDataToPwaMc(int dataentry, int pwamcentry) {
  //  cout << "Distance from Data Event " << dataentry << " to PwaMc Event " << pwamcentry << endl;

  ntdata->GetEntry(dataentry);
  ntpwamc->GetEntry(pwamcentry);

  float distance = sqrt( 
			 pow((datacosThetaKKpipi-pwamccosThetaKKpipi)/2.,2) +
			 pow((dataphiKKpipi-pwamcphiKKpipi)/(2*3.1415),2) +
			 pow((datacosThetaK1pi1pi2-pwamccosThetaK1pi1pi2)/2. ,2) +
			 pow((dataphiK1pi1pi2-pwamcphiK1pi1pi2)/(2*3.1415) ,2) +
			 pow((datacosThetapi1pi2-pwamccosThetapi1pi2)/2. ,2) +
			 pow((dataphipi1pi2-pwamcphipi1pi2)/(2*3.1415) ,2) +
			 pow((datacosThetapi1-pwamccosThetapi1) ,2.) +
			 pow((dataphipi1-pwamcphipi1)/(3.1415) ,2) +
			 pow((datainvmassK1pi1pi2-pwamcinvmassK1pi1pi2)/2.2 ,2) +
			 pow((datainvmasspi1pi2-pwamcinvmasspi1pi2)/2.3 ,2)
			 );
  
  return distance;

}








float GetDistancePwaMcToPwaMc(int pwamcentry1, int pwamcentry2) {
  //  cout << "Distance from PwaMC Event " << dataentry1 << " to PwaMC Event " << dataentry2 << endl;
  
  ntpwamc->GetEntry(pwamcentry1);
  float pwamccosThetaKKpipi_1 = pwamccosThetaKKpipi;
  float pwamcphiKKpipi_1 = pwamcphiKKpipi;
  float pwamccosThetaK1pi1pi2_1 = pwamccosThetaK1pi1pi2;
  float pwamcphiK1pi1pi2_1 = pwamcphiK1pi1pi2;
  float pwamccosThetapi1pi2_1 = pwamccosThetapi1pi2;
  float pwamcphipi1pi2_1 = pwamcphipi1pi2;
  float pwamccosThetapi1_1 = pwamccosThetapi1;
  float pwamcphipi1_1 = pwamcphipi1;
  float pwamcinvmassK1pi1pi2_1 = pwamcinvmassK1pi1pi2;
  float pwamcinvmasspi1pi2_1 = pwamcinvmasspi1pi2;
  ntpwamc->GetEntry(pwamcentry2);
  float pwamccosThetaKKpipi_2 = pwamccosThetaKKpipi;
  float pwamcphiKKpipi_2 = pwamcphiKKpipi;
  float pwamccosThetaK1pi1pi2_2 = pwamccosThetaK1pi1pi2;
  float pwamcphiK1pi1pi2_2 = pwamcphiK1pi1pi2;
  float pwamccosThetapi1pi2_2 = pwamccosThetapi1pi2;
  float pwamcphipi1pi2_2 = pwamcphipi1pi2;
  float pwamccosThetapi1_2 = pwamccosThetapi1;
  float pwamcphipi1_2 = pwamcphipi1;
  float pwamcinvmassK1pi1pi2_2 = pwamcinvmassK1pi1pi2;
  float pwamcinvmasspi1pi2_2 = pwamcinvmasspi1pi2;
  
  float distance = sqrt( 
			 pow((pwamccosThetaKKpipi_1-pwamccosThetaKKpipi_2)/2.,2) +
			 pow((pwamcphiKKpipi_1-pwamcphiKKpipi_2)/(2*3.1415),2) +
			 pow((pwamccosThetaK1pi1pi2_1-pwamccosThetaK1pi1pi2_2)/2. ,2) +
			 pow((pwamcphiK1pi1pi2_1-pwamcphiK1pi1pi2_2)/(2*3.1415) ,2) +
			 pow((pwamccosThetapi1pi2_1-pwamccosThetapi1pi2_2)/2. ,2) +
			 pow((pwamcphipi1pi2_1-pwamcphipi1pi2_2)/(2*3.1415) ,2) +
			 pow((pwamccosThetapi1_1-pwamccosThetapi1_2) ,2) +
			 pow((pwamcphipi1_1-pwamcphipi1_2)/(3.1415) ,2) +
			 pow((pwamcinvmassK1pi1pi2_1-pwamcinvmassK1pi1pi2_2)/2.2 ,2) +
			 pow((pwamcinvmasspi1pi2_1-pwamcinvmasspi1pi2_2)/2.3 ,2)
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



