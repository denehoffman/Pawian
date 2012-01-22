#include <vector>
#include <map>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <TFile>
#include <TH1F>
#include <TNtuple>

TH1F* dataCosThetaKKpipi;
TH1F* mcCosThetaKKpipi;

TH1F* dataPhiKKpipi;
TH1F* mcPhiKKpipi;

TH1F* dataCosThetaKpipi;
TH1F* mcCosThetaKpipi;

TH1F* dataPhiKpipi;
TH1F* mcPhiKpipi;

TH1F* dataCosThetaKpi;
TH1F* mcCosThetaKpi;

TH1F* dataPhiKpi;
TH1F* mcPhiKpi;

TH1F* dataCosThetapi;
TH1F* mcCosThetapi;

TH1F* dataPhipi;
TH1F* mcPhipi;

TH1F* datainvMassKpipi;
TH1F* mcinvMassKpipi;

TH1F* datainvMassKpi;
TH1F* mcinvMassKpi;



void qaCheck(TString fname="bin/gcc-4.1.2/debug/link-static/Psi2STo2K2PiGam.root"){

  using namespace std;
  gROOT->SetStyle("Plain");

  gStyle->SetCanvasColor(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111); 

  std::vector<TH1F*> histVectData;
  std::vector<TH1F*> histVectMc;
  std::vector<TH2F*> histVectData2d;
  std::vector<TH2F*> histVectMc2d;

  TFile* f1 = new TFile(fname,"READ");
  TNtuple* ntdata = (TNtuple*)f1->Get("dataTuple");
  TNtuple* ntmc = (TNtuple*)f1->Get("mcTuple");

  Int_t evtNo = 0;
  float datacosThetaKKpipi;
  float dataphiKKpipi;
  float datacosThetaK1pi1pi2;
  float datacosThetaK2pi1pi2;
  float dataphiK1pi1pi2;
  float dataphiK2pi1pi2;
  float datacosThetaK1pi1;
  float datacosThetaK1pi2;
  float datacosThetaK2pi1;
  float datacosThetaK2pi2;
  float dataphiK1pi1;
  float dataphiK1pi2;
  float dataphiK2pi1;
  float dataphiK2pi2;
  float datacosThetapi1;
  float datacosThetapi2;
  float dataphipi1;
  float dataphipi2;
  float datainvmassK1pi1pi2;
  float datainvmassK2pi1pi2;
  float datainvmassK1pi1;
  float datainvmassK1pi2;  
  float datainvmassK2pi1;  
  float datainvmassK2pi2;
  float dataweight;
  float mccosThetaKKpipi;
  float mcphiKKpipi;
  float mccosThetaK1pi1pi2;
  float mccosThetaK2pi1pi2;
  float mcphiK1pi1pi2;
  float mcphiK2pi1pi2;
  float mccosThetaK1pi1;
  float mccosThetaK1pi2;
  float mccosThetaK2pi1;
  float mccosThetaK2pi2;
  float mcphiK1pi1;
  float mcphiK1pi2;
  float mcphiK2pi1;
  float mcphiK2pi2;
  float mccosThetapi1;
  float mccosThetapi2;
  float mcphipi1;
  float mcphipi2;
  float mcinvmassK1pi1pi2;
  float mcinvmassK2pi1pi2;
  float mcinvmassK1pi1;
  float mcinvmassK1pi2;  
  float mcinvmassK2pi1;  
  float mcinvmassK2pi2;  
  float mcweight;

  ntdata->SetBranchAddress("costhetaKKpipi", &datacosThetaKKpipi);
  ntdata->SetBranchAddress("phiKKpipi", &dataphiKKpipi);

  ntdata->SetBranchAddress("costhetaK1pi1pi2", &datacosThetaK1pi1pi2);
  ntdata->SetBranchAddress("phiK1pi1pi2", &dataphiK1pi1pi2);
  ntdata->SetBranchAddress("costhetaK2pi1pi2", &datacosThetaK2pi1pi2);
  ntdata->SetBranchAddress("phiK2pi1pi2", &dataphiK2pi1pi2);

  ntdata->SetBranchAddress("costhetaK1pi1", &datacosThetaK1pi1);
  ntdata->SetBranchAddress("phiK1pi1", &dataphiK1pi1);
  ntdata->SetBranchAddress("costhetaK1pi2", &datacosThetaK1pi2);
  ntdata->SetBranchAddress("phiK1pi2", &dataphiK1pi2);
  ntdata->SetBranchAddress("costhetaK2pi1", &datacosThetaK2pi1);
  ntdata->SetBranchAddress("phiK2pi1", &dataphiK2pi1);
  ntdata->SetBranchAddress("costhetaK2pi2", &datacosThetaK2pi2);
  ntdata->SetBranchAddress("phiK2pi2", &dataphiK2pi2);

  ntdata->SetBranchAddress("costhetapi1", &datacosThetapi1);
  ntdata->SetBranchAddress("phipi1", &dataphipi1);
  ntdata->SetBranchAddress("costhetapi2", &datacosThetapi2);
  ntdata->SetBranchAddress("phipi2", &dataphipi2);
  
  ntdata->SetBranchAddress("mk1pi1pi2", &datainvmassK1pi1pi2);
  ntdata->SetBranchAddress("mk2pi1pi2", &datainvmassK2pi1pi2);

  ntdata->SetBranchAddress("mk1pi1", &datainvmassK1pi1);
  ntdata->SetBranchAddress("mk1pi2", &datainvmassK1pi2);
  ntdata->SetBranchAddress("mk2pi1", &datainvmassK2pi1);
  ntdata->SetBranchAddress("mk2pi2", &datainvmassK2pi2);

  ntdata->SetBranchAddress("weight", &dataweight);




  ntmc->SetBranchAddress("costhetaKKpipi", &mccosThetaKKpipi);
  ntmc->SetBranchAddress("phiKKpipi", &mcphiKKpipi);

  ntmc->SetBranchAddress("costhetaK1pi1pi2", &mccosThetaK1pi1pi2);
  ntmc->SetBranchAddress("phiK1pi1pi2", &mcphiK1pi1pi2);
  ntmc->SetBranchAddress("costhetaK2pi1pi2", &mccosThetaK2pi1pi2);
  ntmc->SetBranchAddress("phiK2pi1pi2", &mcphiK2pi1pi2);

  ntmc->SetBranchAddress("costhetaK1pi1", &mccosThetaK1pi1);
  ntmc->SetBranchAddress("phiK1pi1", &mcphiK1pi1);
  ntmc->SetBranchAddress("costhetaK1pi2", &mccosThetaK1pi2);
  ntmc->SetBranchAddress("phiK1pi2", &mcphiK1pi2);
  ntmc->SetBranchAddress("costhetaK2pi1", &mccosThetaK2pi1);
  ntmc->SetBranchAddress("phiK2pi1", &mcphiK2pi1);
  ntmc->SetBranchAddress("costhetaK2pi2", &mccosThetaK2pi2);
  ntmc->SetBranchAddress("phiK2pi2", &mcphiK2pi2);
  
  ntmc->SetBranchAddress("costhetapi1", &mccosThetapi1);
  ntmc->SetBranchAddress("phipi1", &mcphipi1);
  ntmc->SetBranchAddress("costhetapi2", &mccosThetapi2);
  ntmc->SetBranchAddress("phipi2", &mcphipi2);
  
  ntmc->SetBranchAddress("mk1pi1pi2", &mcinvmassK1pi1pi2);
  ntmc->SetBranchAddress("mk2pi1pi2", &mcinvmassK2pi1pi2);  
  
  ntmc->SetBranchAddress("mk1pi1", &mcinvmassK1pi1);
  ntmc->SetBranchAddress("mk1pi2", &mcinvmassK1pi2);
  ntmc->SetBranchAddress("mk2pi1", &mcinvmassK2pi1);
  ntmc->SetBranchAddress("mk2pi2", &mcinvmassK2pi2);  

  ntmc->SetBranchAddress("weight", &mcweight);
 
  Int_t numOfEntriesData = ntdata->GetEntries();
  Int_t numOfEntriesMc = ntmc->GetEntries();

  dataCosThetaKKpipi = new TH1F("dataCosThetaKKpipi","dataCosThetaKKpipi",200,-1.,1.);
  dataPhiKKpipi = new TH1F("dataPhiKKpipi","dataPhiKKpipi",200,-4.,4.);
  dataCosThetaKpipi = new TH1F("dataCosThetaKpipi","dataCosThetaKpipi",200,-1.,1.);
  dataPhiKpipi = new TH1F("dataPhiKpipi","dataPhiKpipi",200,-4.,4.);
  dataCosThetaKpi = new TH1F("dataCosThetaKpi","dataCosThetaKpi",200,-1.,1.);
  dataPhiKpi = new TH1F("dataPhiKpi","dataPhiKpi",200,-4.,4.);
  dataCosThetapi = new TH1F("dataCosThetapi","dataCosThetapi",200,-1.,1.);
  dataPhipi = new TH1F("dataPhipi","dataPhipi",200,-4.,4.);
  datainvMassKpipi = new TH1F("datainvMassKpipi", "datainvMassKpipi", 176, 0.8, 3.0);
  datainvMassKpi = new TH1F("datainvMassKpi", "datainvMassKpi", 276, 0.5, 2.8);
   
  mcCosThetaKKpipi = new TH1F("mcCosThetaKKpipi","mcCosThetaKKpipi",200,-1.,1.);
  mcPhiKKpipi = new TH1F("mcPhiKKpipi","mcPhiKKpipi",200,-4.,4.);
  mcCosThetaKpipi = new TH1F("mcCosThetaKpipi","mcCosThetaKpipi",200,-1.,1.);
  mcPhiKpipi = new TH1F("mcPhiKpipi","mcPhiKpipi",200,-4.,4.);
  mcCosThetaKpi = new TH1F("mcCosThetaKpi","mcCosThetaKpi",200,-1.,1.);
  mcPhiKpi = new TH1F("mcPhiKpi","mcPhiKpi",200,-4.,4.);
  mcCosThetapi = new TH1F("mcCosThetapi","mcCosThetapi",200,-1.,1.);
  mcPhipi = new TH1F("mcPhipi","mcPhipi",200,-4.,4.);
  mcinvMassKpipi = new TH1F("mcinvMassKpipi", "mcinvMassKpipi", 176, 0.8, 3.0);
  mcinvMassKpi = new TH1F("mcinvMassKpi", "mcinvMassKpi", 276, 0.5, 2.8);

  histVectData.push_back(dataCosThetaKKpipi);
  histVectData.push_back(dataPhiKKpipi);
  histVectData.push_back(dataCosThetaKpipi);
  histVectData.push_back(dataPhiKpipi);
  histVectData.push_back(dataCosThetaKpi);
  histVectData.push_back(dataPhiKpi);
  histVectData.push_back(dataCosThetapi);
  histVectData.push_back(dataPhipi);
  histVectData.push_back(datainvMassKpipi);
  histVectData.push_back(datainvMassKpi);

  histVectMc.push_back(mcCosThetaKKpipi);
  histVectMc.push_back(mcPhiKKpipi);
  histVectMc.push_back(mcCosThetaKpipi);
  histVectMc.push_back(mcPhiKpipi);
  histVectMc.push_back(mcCosThetaKpi);
  histVectMc.push_back(mcPhiKpi);
  histVectMc.push_back(mcCosThetapi);
  histVectMc.push_back(mcPhipi);
  histVectMc.push_back(mcinvMassKpipi);
  histVectMc.push_back(mcinvMassKpi);

  Int_t maxEvtNo = 0;
  for(Int_t i=0; i<max(numOfEntriesMc,numOfEntriesData); i++)
    {
      if(i<numOfEntriesData){
      	ntdata->GetEntry(i);

	dataCosThetaKKpipi->Fill(datacosThetaKKpipi, dataweight);

        dataPhiKKpipi->Fill(dataphiKKpipi, dataweight);	

        dataCosThetaKpipi->Fill(datacosThetaK1pi1pi2, dataweight);	

        dataPhiKpipi->Fill(dataphiK1pi1pi2, dataweight);	

	dataCosThetaKpi->Fill(datacosThetaK1pi1, dataweight);
	dataCosThetaKpi->Fill(datacosThetaK1pi2, dataweight);  
        
	dataPhiKpi->Fill(dataphiK1pi1, dataweight);                                              	
        dataPhiKpi->Fill(dataphiK1pi2, dataweight);        
        
        dataCosThetapi->Fill(datacosThetapi1, dataweight);                                              	
        dataCosThetapi->Fill(datacosThetapi2, dataweight); 

        dataPhipi->Fill(dataphipi1, dataweight);                                              	
        dataPhipi->Fill(dataphipi2, dataweight);
 
        datainvMassKpipi->Fill(datainvmassK1pi1pi2, dataweight);

        datainvMassKpi->Fill(datainvmassK1pi1, dataweight);
        datainvMassKpi->Fill(datainvmassK1pi2, dataweight);

      }

      ntmc->GetEntry(i);

      mcCosThetaKKpipi->Fill(mccosThetaKKpipi, mcweight);

      mcPhiKKpipi->Fill(mcphiKKpipi, mcweight);

      mcCosThetaKpipi->Fill(mccosThetaK1pi1pi2, mcweight);	

      mcPhiKpipi->Fill(mcphiK1pi1pi2, mcweight);	

      mcCosThetaKpi->Fill(mccosThetaK1pi1, mcweight);                                              	
      mcCosThetaKpi->Fill(mccosThetaK1pi2, mcweight);                                              	

      mcPhiKpi->Fill(mcphiK1pi1, mcweight);                                              	
      mcPhiKpi->Fill(mcphiK1pi2, mcweight);                                              	

      mcCosThetapi->Fill(mccosThetapi1, mcweight);                                              	
      mcCosThetapi->Fill(mccosThetapi2, mcweight); 

      mcPhipi->Fill(mcphipi1, mcweight);                                              	
      mcPhipi->Fill(mcphipi2, mcweight);
      mcinvMassKpipi->Fill(mcinvmassK1pi1pi2, mcweight);

      mcinvMassKpi->Fill(mcinvmassK1pi1, mcweight);
      mcinvMassKpi->Fill(mcinvmassK1pi2, mcweight);

      if(maxEvtNo<evtNo) maxEvtNo=evtNo;
    }


  double scalefactor = double(ntdata->GetEntries())/double(ntmc->GetEntries()); 

  Int_t rebin = 4;

  TCanvas* cmain = new TCanvas("cmain","cmain",1400,600);
  cmain->Divide(5,2);
  for(int i=0; i<histVectData.size(); i++) {
    cmain->cd(i+1);
    histVectData[i]->Rebin(rebin);
    histVectData[i]->SetLineColor(2);
    histVectData[i]->Draw("E");
    histVectMc[i]->Rebin(rebin);
    histVectMc[i]->Scale(scalefactor);
    histVectMc[i]->Draw("same");
    } 


  float chisq = 0;

  float constraints = 98; // 12-01-06_FitParamHyp8_base_K0K2_Hyp5/rerun_correctedampK0K2/Psi2STo2K2PiGamPWA_noK1_1680Hyp7.root
  // float constraints = 94;   // 12-01-06_FitParamHyp8_base_K0K2_Hyp5/rerun_correctedampK0K2/Psi2STo2K2PiGamPWA_noK1_1680Hyp7_noK1_2300Hyp7.root
  float ndf;
  float conflevel = 0;
  int numberofbins = 0;
  
  for(int i=0; i<histVectData.size(); i++) {
    int currentNoOfBins=0;

    for (int j=0; j<(histVectData[i]->GetSize()); j++) {
      float binContentData=histVectData[i]->GetBinContent(j);
      if (binContentData>0.){
	chisq += pow( binContentData - (histVectMc[i]->GetBinContent(j)),2)/binContentData;
	currentNoOfBins++;
      }
      cout << "chi^2 after bin " << j << " of histogram " << histVectData[i]->GetName() << ": " << chisq << endl;
    }
    cout << "chi^2 after " << histVectData[i]->GetName() << ": " << chisq << endl;
    
    numberofbins += currentNoOfBins;
  }   

    ndf = numberofbins - constraints;

    cout << "Total chi^2:\t" << chisq << endl;    
    cout << "Total number of bins:\t" << numberofbins << endl;
    cout << "Number of constraints:\t" << constraints << endl;
    cout << "Degrees of freedom:\t" << ndf << endl;
    cout << "chi^2/ndf:\t" << chisq/ndf << endl;
    cout << "CL:\t" << TMath::Prob(chisq, ndf) << endl;
    cout << endl;

}

 


