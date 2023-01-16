
void drawCrontribs(){
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Divide(2,4);
  c1->cd(1);
  DrawContributions("pawianHists.root","pawianHists.root pawianHistsf2_1270.root pawianHistsa2_1320.root", "Fitpion0Aeta Fitpion0Beta", true, "all f2(1270) a2(1320)");
  c1->cd(2);
  DrawContributions("pawianHists.root","pawianHists.root pawianHistsf2_1270.root pawianHistsa2_1320.root", "Fitpion0Apion0B", true, "all f2(1270) a2(1320");
  c1->cd(3);
  DrawContributions("pawianHists.root","pawianHists.root pawianHistsf2_1270.root pawianHistsa2_1320.root", "FitThetaHeli_eta_Frompion0Apion0Beta", true, "all f2(1270) a2(1320)");
  c1->cd(4);
  DrawContributions("pawianHists.root","pawianHists.root pawianHistsf2_1270.root pawianHistsa2_1320.root", "FitPhiHeli_pion0A_Frompion0Apion0Beta FitPhiHeli_pion0B_Frompion0Apion0Beta", true, "all f2(1270) a2(1320)");
  c1->cd(5);
  DrawContributions("pawianHists.root","pawianHists.root pawianHistsf2_1270.root pawianHistsa2_1320.root", "FitThetaHeli_pion0A_Frompion0Apion0B FitThetaHeli_pion0B_Frompion0Apion0B", true, "all f2(1270) a2(1320)");
  c1->cd(6);
  DrawContributions("pawianHists.root","pawianHists.root pawianHistsf2_1270.root pawianHistsa2_1320.root", "FitPhiHeli_pion0A_Frompion0Apion0B FitPhiHeli_pion0B_Frompion0Apion0B", true, "all f2(1270) a2(1320)");
  c1->cd(7);
  DrawContributions("pawianHists.root","pawianHists.root pawianHistsf2_1270.root pawianHistsa2_1320.root", "FitThetaHeli_eta_Frompion0Aeta FitThetaHeli_eta_Frompion0Beta", true, "all f2(1270) a2(1320)");
  c1->cd(8);
  DrawContributions("pawianHists.root","pawianHists.root pawianHistsf2_1270.root pawianHistsa2_1320.root", "FitPhiHeli_eta_Frompion0Aeta FitPhiHeli_eta_Frompion0Beta", true, "all f2(1270) a2(1320)");
  c1->SaveAs("FitResultHists.jpg");  
}
