TFile *f;
int ipad;
TCanvas* csingle;

void plotQA(TString fName="./JpsiGamPipPimPi0KK.root"){
  
    gStyle->SetPadBorderMode(1);
  gStyle->SetPadBorderSize(1);
  gStyle->SetPadColor(0);
  gStyle->SetDrawBorder(0);
  gStyle->SetCanvasBorderMode(1);
  gStyle->SetCanvasBorderSize(1);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);

  gStyle->SetTitleColor(1);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleY(0.975);
  
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.8,"y");
  gStyle->SetTitleOffset(1.7,"z");
  
  //gStyle->SetLabelOffset(.02,"x");
  gStyle->SetLabelOffset(.025,"Y");
  //gStyle->SetLabelOffset(.02,"z");
    
  gStyle->SetPadTopMargin(0.1); 
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.18);


  gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);


  f = new TFile(fName, "READ");
  if(f==0){
    cout << "Could not open file " << fName << endl;
    return;
  }
  
  csingle=new TCanvas("cs","s1",450,450);
  

  ipad=1;
  TCanvas* c1 = new TCanvas("c1","c1",450,600);
  c1->Divide(2,3);
  TCanvas* c2 = new TCanvas("c2","c2",450,600);
  c2->Divide(2,3);
  
  ipad=1;
  
  c1->cd(ipad);
  plotHist("_OmegaPhiMass", "m(#omega#phi) [GeV/c^{2}]", "Entries / %f GeV");
  
  c1->cd(ipad);
  plotHist("_costGamCm", "cos(#theta^{*}_{#gamma})");
  
  c1->cd(ipad);
  plotHist("_costPip_PipPimPi0Heli", "cos(#theta_{#pi^{+}})" );
  
  c1->cd(ipad);
  plotHist("_costKp_KpKmHeli", "cos(#theta_{K^{+}})");
  
  c1->cd(ipad);
  plotHist("_phiPip_PipPimPi0Heli", "#phi_{#pi^{+}}/#pi");
  
  c1->cd(ipad);
  plotHist("_phiKp_KpKmHeli", "#phi_{K^{+}}/#pi");
  
  

  ipad=1;
  c2->cd(ipad);
  plotHist("_costPhi_OmegaPhiHeli", "cos(#theta_{K^{+}K^{-}})");
  
  c2->cd(ipad);
  plotHist("_phiPhi_OmegaPhiHeli", "#phi_{K^{+}K^{-}}");
  
  c2->cd(ipad);
  plotHist("_chi", "#chi_{#phi#phi} [deg]");
  
  c2->cd(ipad);
  plotHist("lambda","#lambda_{#omega}");

  
  c1->Print("qaPlots.ps(");
  c2->Print("qaPlots.ps)");
  
}



void plotHist(TString base, TString xtitle="", TString ytitle=""){
  
  TH1F* hd = (TH1F*) f->Get(base+"DataHist");
  TH1F* hf = (TH1F*) f->Get(base+"FittedHist");
  
  hd->SetMinimum(0);
  hf->SetMinimum(0);
  
  hd->SetXTitle(xtitle);
  char  tmp[1000];
  sprintf(tmp,"Entries / %.3f", hd->GetBinWidth(1) );
  if(ytitle!=""){
    hd->SetYTitle(ytitle);
  }
  hd->SetYTitle(tmp);
  
  
  
  hf->SetLineWidth(2);
  hf->SetLineColor(kRed);
  
  hd->Draw("E");
  hf->Draw("same");
  
  ipad++;
  
  csingle->cd();
  hd->Draw("E");
  hf->Draw("same");
  csingle->SaveAs(base+".eps");
  return;
}
