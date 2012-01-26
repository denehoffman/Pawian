TFile *f;
int ipad;
void plotQA(TString fName="./JpsiGamKsKlKK.root"){
  
  f = new TFile(fName, "READ");
  if(f==0){
    cout << "Could not open file " << fName << endl;
    return;
  }
  
  ipad=1;
  TCanvas* c1 = new TCanvas("c1","c1",450,600);
  c1->Divide(2,3);
  TCanvas* c2 = new TCanvas("c2","c2",450,600);
  c2->Divide(2,3);
  
  ipad=1;
  
  c1->cd(ipad);
  plotHist("_PhiPhiMass");
  
  c1->cd(ipad);
  plotHist("_costGamCm");
  
  c1->cd(ipad);
  plotHist("_costKs_KsKlHeli");
  
  c1->cd(ipad);
  plotHist("_costKp_KpKmHeli");
  
  c1->cd(ipad);
  plotHist("_phiKs_KsKlHeli");
  
  c1->cd(ipad);
  plotHist("_phiKp_KpKmHeli");
  
  

  ipad=1;
  c2->cd(ipad);
  plotHist("_costPhi_KpKm");
  
  c2->cd(ipad);
  plotHist("_phiPhi_KpKm");
  
  c2->cd(ipad);
  plotHist("_chi");
  
  c1->Print("qaPlots.pdf(");
  c2->Print("qaPlots.pdf)");
  
}



void plotHist(TString base){
  TH1F* hd = (TH1F*) f->Get(base+"DataHist");
  TH1F* hf = (TH1F*) f->Get(base+"FittedHist");
  
  hd->SetMinimum(0);
  hf->SetMinimum(0);
  
  hf->SetLineWidth(2);
  hf->SetLineColor(kRed);
  
  hd->Draw("E");
  hf->Draw("same");
  
  ipad++;
  return;
}
